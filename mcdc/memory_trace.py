"""Opt-in memory diagnostics for the Slurm step on MPI rank zero's node.

The sampler runs in a separate, standard-library-only interpreter so samples
continue during compiled transport and blocking MPI collectives. Records are
appended as JSON lines. Linux counters use bytes; unavailable values are null.
"""

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import select
import socket
import subprocess
import sys
import time
import uuid

try:
    import fcntl
except ImportError:
    fcntl = None


_trace = None


def _read_number(path):
    try:
        return int(path.read_text().strip())
    except (OSError, ValueError):
        return None


def _unescape_mount(value):
    return re.sub(r"\\([0-7]{3})", lambda match: chr(int(match[1], 8)), value)


def _find_step_cgroup(proc_root=Path("/proc")):
    """Resolve an accessible Slurm step ancestor in a memory cgroup mount."""
    try:
        memberships = (proc_root / "self/cgroup").read_text().splitlines()
        mount_lines = (proc_root / "self/mountinfo").read_text().splitlines()
    except OSError:
        return None

    for membership in memberships:
        fields = membership.split(":", 2)
        if len(fields) != 3:
            continue
        _, controllers, group = fields
        if "memory" in controllers.split(","):
            version = 1
        elif controllers == "":
            version = 2
        else:
            continue
        for line in mount_lines:
            try:
                before, after = line.split(" - ", 1)
                mount_fields = before.split()
                fs_type, _, options = after.split()[:3]
                root = Path(_unescape_mount(mount_fields[3]))
                mount = Path(_unescape_mount(mount_fields[4]))
            except (ValueError, IndexError):
                continue
            if version == 1 and (
                fs_type != "cgroup" or "memory" not in options.split(",")
            ):
                continue
            if version == 2 and fs_type != "cgroup2":
                continue
            try:
                relative = Path(group).relative_to(root)
            except ValueError:
                continue
            current = mount / relative
            for candidate in (current, *current.parents):
                name = root.name if candidate == mount else candidate.name
                if re.fullmatch(r"step[_-][^/]+", name):
                    counter = (
                        "memory.current" if version == 2 else "memory.usage_in_bytes"
                    )
                    if (candidate / counter).is_file():
                        return {
                            "path": str(candidate),
                            "mount": str(mount),
                            "version": version,
                        }
                if candidate == mount:
                    break
    return None


def _cgroup_memory(group):
    result = {
        "cgroup_current_bytes": None,
        "cgroup_peak_bytes": None,
        "cgroup_limit_bytes": None,
        "cgroup_limit_path": None,
    }
    if group is None:
        return result
    path = Path(group["path"])
    mount = Path(group["mount"])
    if group["version"] == 2:
        current_name, peak_name, limit_name = (
            "memory.current",
            "memory.peak",
            "memory.max",
        )
    else:
        current_name, peak_name, limit_name = (
            "memory.usage_in_bytes",
            "memory.max_usage_in_bytes",
            "memory.limit_in_bytes",
        )
    result["cgroup_current_bytes"] = _read_number(path / current_name)
    result["cgroup_peak_bytes"] = _read_number(path / peak_name)
    # Slurm often limits the parent job while leaving the step unlimited.
    # Report the tightest visible ancestor limit with its source path.
    for parent in (path, *path.parents):
        value = _read_number(parent / limit_name)
        if value is not None and 0 <= value < 2**60:
            previous = result["cgroup_limit_bytes"]
            if previous is None or value < previous:
                result["cgroup_limit_bytes"] = value
                result["cgroup_limit_path"] = str(parent / limit_name)
        if parent == mount:
            break
    return result


def _rss(pid, proc_root=Path("/proc")):
    if pid is None:
        return None
    try:
        for line in (proc_root / str(pid) / "status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) * 1024
    except (OSError, ValueError, IndexError):
        pass
    return None


def _write_record(fd, context, event, **fields):
    record = {
        "time": datetime.now(timezone.utc).isoformat(),
        "monotonic_ns": time.monotonic_ns(),
        "run_id": context["run_id"],
        "event": event,
        **_cgroup_memory(context["cgroup"]),
        "rss_rank0_bytes": _rss(context["rank0_pid"]),
        "rss_rank1_bytes": _rss(context["worker_pid"]),
        **fields,
    }
    payload = (json.dumps(record, separators=(",", ":")) + "\n").encode()
    # Parent and sampler use separate open-file descriptions. Locking keeps
    # each JSON line intact even on filesystems without atomic append writes.
    fcntl.flock(fd, fcntl.LOCK_EX)
    try:
        while payload:
            written = os.write(fd, payload)
            if written == 0:
                raise OSError("Memory trace write returned zero bytes")
            payload = payload[written:]
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)


def _warning(error):
    print(f"MC/DC memory trace unavailable: {error}", file=sys.stderr, flush=True)


def start(path, worker_pid=None):
    """Start a rank-zero trace, returning whether markers can be recorded.

    Call on rank zero. ``worker_pid`` is a rank-one process on the same host.
    The file is appended so each invocation retains a distinct ``run_id``.
    Sampling or filesystem failures leave the simulation running.
    """
    global _trace
    stop()
    fd = None
    try:
        if fcntl is None:
            raise OSError("POSIX file locking is required")
        path = Path(path).expanduser().resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o600)
        context = {
            "run_id": uuid.uuid4().hex,
            "rank0_pid": os.getpid(),
            "worker_pid": worker_pid,
            "cgroup": _find_step_cgroup(),
        }
        _write_record(
            fd,
            context,
            "trace_start",
            hostname=socket.gethostname(),
            rank0_pid=context["rank0_pid"],
            rank1_pid=worker_pid,
            slurm_job_id=os.environ.get("SLURM_JOB_ID"),
            slurm_step_id=os.environ.get("SLURM_STEP_ID"),
            cgroup=context["cgroup"],
            sample_interval_seconds=1.0,
            scope="rank-zero-node Slurm step",
        )
        _trace = {"fd": fd, "context": context, "process": None}
    except (OSError, ValueError, TypeError) as error:
        if fd is not None:
            os.close(fd)
        _warning(error)
        return False
    print(f"MC/DC memory trace: {path}", flush=True)
    if context["cgroup"] is None:
        _warning("Slurm step cgroup counters are unavailable on this host")
    try:
        # Run the file directly rather than importing mcdc in the sampler.
        _trace["process"] = subprocess.Popen(
            [sys.executable, "-I", __file__, str(path), json.dumps(context)],
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, ValueError) as error:
        mark("sampler_unavailable", tally={"error": str(error)})
        _warning(error)
    return True


def mark(event, batch=None, tally=None):
    """Write a timestamped operation marker when a trace is active."""
    global _trace
    if _trace is None:
        return
    try:
        process = _trace["process"]
        if process is not None and process.poll() is not None:
            _trace["process"] = None
            process.stdin.close()
            _write_record(
                _trace["fd"],
                _trace["context"],
                "sampler_exited",
                returncode=process.returncode,
            )
            _warning(f"sampler exited with status {process.returncode}")
        _write_record(_trace["fd"], _trace["context"], event, batch=batch, tally=tally)
    except (OSError, ValueError, TypeError) as error:
        _warning(error)
        stop()


def stop():
    """Finish the current trace and shut down its independent sampler."""
    global _trace
    if _trace is None:
        return
    trace, _trace = _trace, None
    process = trace["process"]
    if process is not None:
        try:
            process.stdin.close()
            process.wait(timeout=2)
        except (OSError, subprocess.TimeoutExpired):
            try:
                process.kill()
                process.wait(timeout=1)
            except (OSError, subprocess.TimeoutExpired):
                pass
    try:
        _write_record(trace["fd"], trace["context"], "trace_stop")
    except (OSError, ValueError, TypeError):
        pass
    finally:
        try:
            os.close(trace["fd"])
        except OSError:
            pass


def _sample(path, context):
    """Sample until rank zero closes its pipe or dies, including SIGKILL."""
    fd = os.open(path, os.O_WRONLY | os.O_APPEND)
    try:
        while True:
            _write_record(fd, context, "sample")
            # EOF arrives automatically if rank zero is OOM-killed. The child
            # never imports MPI and never participates in a collective.
            readable, _, _ = select.select([sys.stdin], [], [], 1.0)
            if readable:
                break
        _write_record(fd, context, "sampler_stop")
    finally:
        os.close(fd)


if __name__ == "__main__":
    try:
        _sample(sys.argv[1], json.loads(sys.argv[2]))
    except (OSError, ValueError, IndexError):
        sys.exit(1)
