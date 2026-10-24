"""Run solver.py against the disclosed practice case and print metrics with known truth."""
import json
import sys
from pathlib import Path

import harness
import lab

HERE = Path(__file__).parent


def main():
    case = json.loads((HERE / "practice_case.json").read_text())
    L = lab.Lab(case)
    out = harness.run_session(L, str(HERE / "solver.py"), str(HERE))
    if out["error"]:
        print("session failed:", out["error"])
        return 1
    m = lab.evaluate(L, out["result"])
    print("budget spent: %.2f of %.0f; donors sequenced: %s; batches used: %s" % (m["spent"], lab.BUDGET, m["donors_sequenced"], m["batches_used"]))
    print("true effects called with the right sign: %d of %d; false calls: %d; true effects inside intervals: %d; null interval coverage: %.2f; tests missing: %d"
          % (m["true_called_correct_sign"], m["true_effects"], m["false_calls"], m["true_covered"], m["null_coverage"], m["missing_tests"]))
    ok, why = lab.check_rules(m)
    print("meets the published rules on this case:", ok, why)
    return 0


if __name__ == "__main__":
    sys.exit(main())
