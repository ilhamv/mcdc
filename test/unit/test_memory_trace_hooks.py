"""Exercise diagnostic lifetime without launching a sampler or MPI workers."""

from types import SimpleNamespace
from unittest.mock import Mock

import pytest
from mpi4py import MPI

import mcdc.config as config
from mcdc.main import trace_run_memory
from mcdc import memory_trace


def test_trace_stops_after_transport_exception(monkeypatch):
    monkeypatch.setattr(config.args, "memory_trace", "memory.jsonl")
    comm = SimpleNamespace(rank=0, size=2, bcast=Mock(return_value=("host", 123)))
    monkeypatch.setattr(MPI, "COMM_WORLD", comm)
    monkeypatch.setattr(MPI, "Get_processor_name", lambda: "host")
    start = Mock()
    stop = Mock()
    monkeypatch.setattr(memory_trace, "start", start)
    monkeypatch.setattr(memory_trace, "stop", stop)

    with pytest.raises(RuntimeError, match="transport failed"):
        with trace_run_memory():
            raise RuntimeError("transport failed")

    start.assert_called_once_with("memory.jsonl", worker_pid=123)
    stop.assert_called_once_with()


def test_disabled_trace_has_no_mpi_or_sampler_activity(monkeypatch):
    monkeypatch.setattr(config.args, "memory_trace", None)
    monkeypatch.setattr(MPI, "COMM_WORLD", None)
    start = Mock()
    stop = Mock()
    monkeypatch.setattr(memory_trace, "start", start)
    monkeypatch.setattr(memory_trace, "stop", stop)

    with trace_run_memory():
        pass

    start.assert_not_called()
    stop.assert_not_called()


def test_other_node_worker_is_not_read_through_local_proc(monkeypatch):
    monkeypatch.setattr(config.args, "memory_trace", "memory.jsonl")
    comm = SimpleNamespace(rank=0, size=2, bcast=Mock(return_value=("other", 123)))
    monkeypatch.setattr(MPI, "COMM_WORLD", comm)
    monkeypatch.setattr(MPI, "Get_processor_name", lambda: "host")
    start = Mock()
    monkeypatch.setattr(memory_trace, "start", start)
    monkeypatch.setattr(memory_trace, "stop", Mock())

    with trace_run_memory():
        pass

    start.assert_called_once_with("memory.jsonl", worker_pid=None)
