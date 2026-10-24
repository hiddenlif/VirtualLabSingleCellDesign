"""Serves a simulated lab to a solver process over newline-delimited JSON.

The solver runs as a child process and sees only the messages it is sent. The simulator
state, including any hidden truth, stays in this parent process.
"""

import json
import os
import queue
import signal
import subprocess
import sys
import threading
import time


class LabError(Exception):
    """An action the instrument refuses; reported to the solver as ok=false."""


MAX_LINE = 1 << 20  # longest accepted protocol line, in characters
OVERSIZE = object()


def _reader(stream, q):
    while True:
        line = stream.readline(MAX_LINE + 1)
        if not line:
            break
        if len(line) > MAX_LINE:
            q.put(OVERSIZE)
            return
        q.put(line)
    q.put(None)


def _drain(stream, sink):
    while True:
        line = stream.readline(4096)  # bounded reads; long stderr lines are truncated
        if not line:
            break
        if len(sink) < 50:
            sink.append(line.rstrip()[:300])


def _kill_uid(uid):
    """SIGKILL every remaining process owned by uid (grader mode only; Linux /proc)."""
    if not os.path.isdir("/proc"):
        return
    for name in os.listdir("/proc"):
        if not name.isdigit():
            continue
        try:
            with open(f"/proc/{name}/status") as fh:
                owner = next(
                    int(line.split()[1]) for line in fh if line.startswith("Uid:")
                )
            if owner == uid:
                os.kill(int(name), signal.SIGKILL)
        except Exception:
            pass


def run_session(lab, solver_path, cwd, wall_seconds=900, max_messages=5000, uid=None):
    """Run one solver against one lab instance.

    Returns {"result": submitted dict or None, "error": str or None, "messages": int,
    "trace": list of (action, ok) pairs}.
    """
    demote = None
    privileged = uid is not None and hasattr(os, "geteuid") and os.geteuid() == 0
    if privileged:

        def demote():
            os.setgroups([])  # drop root's supplementary groups before changing ids
            os.setgid(uid)
            os.setuid(uid)

    env = {
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
        "PYTHONHASHSEED": "0",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
    proc = subprocess.Popen(
        [sys.executable, "-B", "-s", "-u", solver_path],
        cwd=cwd,
        env=env,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
        preexec_fn=demote,
        start_new_session=True,
    )
    out_q = queue.Queue()
    err_lines = []
    threading.Thread(target=_reader, args=(proc.stdout, out_q), daemon=True).start()
    threading.Thread(target=_drain, args=(proc.stderr, err_lines), daemon=True).start()
    deadline = time.time() + wall_seconds
    trace = []
    outcome = {"result": None, "error": None, "messages": 0, "trace": trace}

    def send(obj):
        proc.stdin.write(json.dumps(obj) + "\n")
        proc.stdin.flush()

    try:
        send({"type": "hello", "protocol": 1, "config": lab.public_config()})
        while True:
            remaining = deadline - time.time()
            if remaining <= 0:
                outcome["error"] = "wall-clock limit reached"
                break
            try:
                line = out_q.get(timeout=min(remaining, 5.0))
            except queue.Empty:
                continue
            if line is OVERSIZE:
                outcome["error"] = f"protocol line longer than {MAX_LINE} characters"
                break
            if line is None:
                outcome["error"] = "solver exited without submitting" + (
                    ": " + " | ".join(err_lines[-3:]) if err_lines else ""
                )
                break
            outcome["messages"] += 1
            if outcome["messages"] > max_messages:
                outcome["error"] = "message limit reached"
                break
            try:
                msg = json.loads(line)
                action = msg.pop("action")
            except Exception:
                outcome["error"] = "malformed message"
                break
            if action == "submit":
                outcome["result"] = msg.get("result")
                trace.append((action, True))
                break
            try:
                reply = {"ok": True, "result": lab.handle(action, msg)}
                trace.append((action, True))
            except LabError as e:
                reply = {"ok": False, "error": str(e)}
                trace.append((action, False))
            except Exception as e:  # a malformed request must not crash the harness
                reply = {"ok": False, "error": f"invalid request: {type(e).__name__}"}
                trace.append((action, False))
            send(reply)
    except BrokenPipeError:
        outcome["error"] = outcome["error"] or "solver closed the pipe"
    finally:
        try:
            proc.stdin.close()
        except Exception:
            pass
        try:
            os.killpg(
                proc.pid, signal.SIGKILL
            )  # the solver and any children it started
        except Exception:
            pass
        proc.kill()
        try:
            proc.wait(timeout=5)
        except Exception:
            pass
        if privileged:
            _kill_uid(uid)  # processes that left the session group
    return outcome
