"""Memory accounting must identify the Slurm step and its inherited limits."""

import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time

import pytest

from mcdc import memory_trace


@pytest.fixture
def cgroup_files(tmp_path):
    def create(version, mount_root="/", mount_name="cgroup", group=None):
        group = group or "/slurm/uid_1/job_42/step_0/task_0"
        proc = tmp_path / "proc"
        (proc / "self").mkdir(parents=True, exist_ok=True)
        mount = tmp_path / mount_name
        current = mount / Path(group).relative_to(mount_root)
        current.mkdir(parents=True, exist_ok=True)
        controllers = "memory" if version == 1 else ""
        fs_type = "cgroup" if version == 1 else "cgroup2"
        super_options = "rw,memory" if version == 1 else "rw"
        mount_text = str(mount).replace(" ", r"\040")
        (proc / "self/cgroup").write_text(f"0:{controllers}:{group}\n")
        (proc / "self/mountinfo").write_text(
            f"30 25 0:28 {mount_root} {mount_text} rw - "
            f"{fs_type} cgroup {super_options}\n"
        )
        return proc, mount, current

    return create


def test_v2_reads_step_counters_and_inherited_job_limit(cgroup_files):
    proc, mount, task = cgroup_files(2)
    step, job = task.parent, task.parent.parent
    (step / "memory.current").write_text("123\n")
    (step / "memory.peak").write_text("456\n")
    (step / "memory.max").write_text("max\n")
    (job / "memory.max").write_text("1000\n")
    (mount / "memory.max").write_text("2000\n")
    # A task's smaller accounting scope must never replace the step counters.
    (task / "memory.current").write_text("7\n")
    group = memory_trace._find_step_cgroup(proc)
    assert group["path"] == str(step)
    assert memory_trace._cgroup_memory(group) == {
        "cgroup_current_bytes": 123,
        "cgroup_peak_bytes": 456,
        "cgroup_limit_bytes": 1000,
        "cgroup_limit_path": str(job / "memory.max"),
    }


def test_v1_handles_mount_roots_escapes_and_unlimited_sentinel(cgroup_files):
    proc, _, task = cgroup_files(1, mount_root="/slurm", mount_name="memory group")
    step, job = task.parent, task.parent.parent
    (step / "memory.usage_in_bytes").write_text("400\n")
    (step / "memory.max_usage_in_bytes").write_text("700\n")
    (step / "memory.limit_in_bytes").write_text(str(2**63 - 4096))
    (job / "memory.limit_in_bytes").write_text("2000\n")
    group = memory_trace._find_step_cgroup(proc)
    assert group["path"] == str(step)
    assert group["version"] == 1
    values = memory_trace._cgroup_memory(group)
    assert values["cgroup_current_bytes"] == 400
    assert values["cgroup_peak_bytes"] == 700
    assert values["cgroup_limit_bytes"] == 2000


def test_step_can_be_the_memory_mount_root(cgroup_files):
    root = "/slurm/uid_1/job_42/step_0"
    proc, mount, _ = cgroup_files(2, mount_root=root)
    (mount / "memory.current").write_text("100\n")
    assert memory_trace._find_step_cgroup(proc)["path"] == str(mount)


def test_does_not_substitute_job_or_machine_counters(cgroup_files):
    proc, mount, task = cgroup_files(2)
    (task.parent.parent / "memory.current").write_text("1000\n")
    (mount / "memory.current").write_text("10000\n")
    assert memory_trace._find_step_cgroup(proc) is None


def test_missing_peak_and_unlimited_ancestors_are_null(cgroup_files):
    proc, mount, task = cgroup_files(2)
    step = task.parent
    (step / "memory.current").write_text("123\n")
    (step / "memory.max").write_text("max\n")
    (mount / "memory.max").write_text("max\n")
    values = memory_trace._cgroup_memory(memory_trace._find_step_cgroup(proc))
    assert values["cgroup_peak_bytes"] is None
    assert values["cgroup_limit_bytes"] is None
    assert values["cgroup_limit_path"] is None


def test_proc_rss_is_converted_to_bytes(tmp_path):
    process = tmp_path / "42"
    process.mkdir()
    (process / "status").write_text("Name:\tpython\nVmRSS:\t1234 kB\n")
    assert memory_trace._rss(42, tmp_path) == 1234 * 1024
    assert memory_trace._rss(43, tmp_path) is None


def test_import_and_disabled_mark_work_without_fcntl(monkeypatch):
    monkeypatch.setitem(sys.modules, "fcntl", None)
    spec = importlib.util.spec_from_file_location(
        "memory_trace_without_fcntl", memory_trace.__file__
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.mark("disabled")
    assert module.start("unused-memory-trace.jsonl") is False
    assert module._trace is None


def test_unwritable_destination_fails_without_starting_sampler(tmp_path):
    assert memory_trace.start(tmp_path) is False
    assert memory_trace._trace is None


@pytest.mark.skipif(memory_trace.fcntl is None, reason="POSIX sampler")
def test_write_failure_stops_trace_without_raising(tmp_path, monkeypatch):
    memory_trace.start(tmp_path / "trace.jsonl")

    def fail(*args, **kwargs):
        raise OSError("simulated full disk")

    monkeypatch.setattr(memory_trace, "_write_record", fail)
    memory_trace.mark("batch_start", batch=1)
    assert memory_trace._trace is None
    memory_trace.stop()


@pytest.mark.skipif(memory_trace.fcntl is None, reason="POSIX sampler")
def test_markers_append_sessions_and_accept_tally_metadata(tmp_path):
    path = tmp_path / "trace.jsonl"
    run_ids = set()
    for _ in range(2):
        try:
            assert memory_trace.start(path)
            memory_trace.mark("before_reduce", batch=1, tally={"id": 0, "bins": 100})
        finally:
            memory_trace.stop()
    records = [json.loads(line) for line in path.read_text().splitlines()]
    for record in records:
        run_ids.add(record["run_id"])
    assert len(run_ids) == 2
    markers = [r for r in records if r["event"] == "before_reduce"]
    assert len(markers) == 2
    assert all(r["tally"] == {"id": 0, "bins": 100} for r in markers)
    assert sum(r["event"] == "trace_stop" for r in records) == 2


@pytest.mark.skipif(memory_trace.fcntl is None, reason="POSIX sampler")
def test_sampler_exits_when_rank_zero_is_killed(tmp_path):
    path = tmp_path / "trace.jsonl"
    driver = """
import importlib.util, sys
spec = importlib.util.spec_from_file_location('memory_trace', sys.argv[1])
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
assert module.start(sys.argv[2])
sys.stdin.read()
"""
    process = subprocess.Popen(
        [sys.executable, "-I", "-c", driver, memory_trace.__file__, str(path)],
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    try:
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            if path.exists() and '"event":"sample"' in path.read_text():
                break
            if process.poll() is not None:
                pytest.fail(process.stderr.read().decode())
            time.sleep(0.05)
        else:
            pytest.fail("Sampler did not start")
        process.kill()
        process.wait(timeout=5)
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            if '"event":"sampler_stop"' in path.read_text():
                break
            time.sleep(0.05)
        else:
            pytest.fail("Sampler did not observe parent-pipe EOF")
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        process.stdin.close()
        process.stderr.close()
