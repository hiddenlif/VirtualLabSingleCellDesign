"""Starter solver. Replace plan() with your own study design and analysis.

The instrument is driven through `inst.call(action, **kwargs)`; see SPEC.md.
"""

from labclient import Instrument


def plan(inst):
    """Return the result dict described in SPEC.md ("Submission")."""
    raise NotImplementedError


if __name__ == "__main__":
    instrument = Instrument()
    instrument.submit(plan(instrument))
