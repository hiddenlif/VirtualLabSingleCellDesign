# Virtual Lab Task: Single-Cell Study Design

Spend a fixed budget on a single-cell RNA follow-up study: recruit donors into batches, choose cells per library, add libraries or resequence, then analyse the returned count matrices to find which genes respond to a treatment in which cell type. Choices change the simulated data; donors, batches, cells and reads are separate resources with different information value.

**Task type:** Interactive simulated experimentation. **Version:** 1.0.0.

All donors, libraries and counts are synthetic. The model is a simplified 40-gene negative-binomial panel with donor and batch effects; it is not a validated power calculator or a realistic transcriptome.

## Contents

- `TASK.md`: task, deliverable and acceptance rules.
- `app/SPEC.md`: protocol and model. `app/lab.py`: simulator, ledger and scoring. `app/harness.py`, `app/labclient.py`: the JSON-lines protocol.
- `app/solver.py`: starter to replace. `app/practice.py`, `app/practice_case.json`: practice case with known truth.
- `SOURCES.md`: source-to-mechanism table, invented fixtures and modeling limitations.
- `Dockerfile`: dependency environment.

## Use

```sh
docker build -t virtual-lab-single-cell-design .
docker run --rm virtual-lab-single-cell-design python3 /app/practice.py
```

The starter raises `NotImplementedError`; the practice runner reports metrics once `solver.py` submits a result.

## Execution evidence

In a local Docker run of the grader, a reference solution was accepted on 8 of 8 held-out cases (reward 1) and the unchanged starter scored 0. These are local simulator results. No blind solver trial or external review has been run, and nothing here is validated against laboratory data. A fixed acquisition plan is a valid solution, so the task does not require adaptivity.

## Reference and verifier

The reference solution, held-out cases, grader, comparison routes and correctness checks are kept in a separate private repository. They are not part of this repository.
