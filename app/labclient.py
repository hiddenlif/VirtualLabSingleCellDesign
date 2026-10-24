"""Client side of the instrument protocol. Import this from solver.py."""
import base64
import json
import sys

try:
    import numpy as np
except ImportError:  # numpy is only needed to decode arrays
    np = None


class LabFault(Exception):
    """The instrument refused an action."""


def decode_array(obj):
    """Decode {"shape": [...], "dtype": "float32", "b64": "..."} into a numpy array."""
    raw = base64.b64decode(obj["b64"])
    return np.frombuffer(raw, dtype=obj["dtype"]).reshape(obj["shape"]).astype("float64")


class Instrument:
    def __init__(self):
        hello = json.loads(sys.stdin.readline())
        self.config = hello["config"]

    def call(self, action, **kwargs):
        sys.stdout.write(json.dumps(dict(kwargs, action=action)) + "\n")
        sys.stdout.flush()
        line = sys.stdin.readline()
        if not line:
            raise LabFault("instrument closed the connection")
        reply = json.loads(line)
        if not reply["ok"]:
            raise LabFault(reply["error"])
        return reply["result"]

    def submit(self, result):
        sys.stdout.write(json.dumps({"action": "submit", "result": result}) + "\n")
        sys.stdout.flush()
