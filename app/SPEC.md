# Single-cell study design virtual lab: specification

All data are synthetic. The graded cases are other draws from `lab.make_case`.

## Protocol

JSON-lines exchange as in `labclient.py`: a `hello` line carries the config; the solver sends
`{"action": ..., ...}` lines, reads `{"ok": true, "result": ...}` or `{"ok": false, "error": ...}`,
and ends with `{"action": "submit", "result": {...}}`.

### Actions

- `status`: budget left, donors available, libraries per batch and each library's donor, condition, batch, cells and depth. Changes nothing.
- `counts` with `library`: the library's current counts (`counts`: array object, cells by 40 genes, `labclient.decode_array`), `condition`, `batch`, `donor`, and QC: median UMI per cell and a noisy `ambient_fraction_estimate`. Free. Refused for an unsequenced library.
- `recruit` with `condition` (`case` or `control`), `batch` (1 to 4), `cells` (200 to 3000): a new donor and its first library. Cost 6 + cells/500. There are 12 new donors of each condition. A batch holds at most 4 libraries.
- `new_library` with `donor`, `batch`, `cells`: another library (fresh cells, same donor). Cost 3 + cells/500.
- `sequence` with `library`, `reads_per_cell` (1000 to 40000): adds reads; cumulative depth is limited to 60000. Cost cells x reads / 5e6. Sequencing again observes more of the same molecules in the same cells.

## Model

A 40-gene panel; G01-G05, G06-G10 and G11-G15 mark type_1, type_2 and type_3 (own-type marker mean 14, other types 4);
G16-G40 are test genes. For a cell of type t from donor d in batch b:

- mean = base[t,g] x exp(donor[d,t,g]) x exp(batch[b,g]) x capture[b] x size, with donor and batch effects log-normal (standard deviations 0.15 and 0.25), size log-normal (0.25)
- cases: mean multiplied by 2^lfc[t,g]; lfc is non-zero for 8 test genes, all in type_3 (magnitudes 0.7 to 1.2)
- true molecules ~ negative binomial (dispersion 0.15), plus ambient molecules whose fraction is 0 to 0.12 in a normal library and 0.3 to 0.5 in a low-quality one
- sequencing to r reads per cell observes each molecule with probability 1 - exp(-r / 8000)
- donors differ in cell-type proportions (about 62%, 30%, 8%)

Pilot: libraries L1 to L4 (donors P1 to P4), 500 cells each, sequenced to 4000 reads per cell.

## Submission

```
{"effects": {"type_1": {"G16": {"log2fc": x, "ci": [lo, hi], "significant": bool}, ...}, "type_2": {...}, "type_3": {...}},
 "ledger": {"spent": x, "libraries": ["L1", ...], ...},
 "follow_up": {"next_step": "recruit_donors" | "new_libraries" | "resequence" | "no_further_acquisition",
               "cell_type": "type_1" | "type_2" | "type_3" | null, "rationale": "..."}}
```

`ledger.spent` must equal the budget used within 0.01 and `ledger.libraries` must list exactly the sequenced libraries (pilot included). Other ledger fields are free. The follow-up is checked for form only.

## Scoring quantities (`lab.evaluate`, `lab.check_rules`)

True effects called with the right sign, false calls among non-responding pairs, true effects inside their intervals, coverage of 0 among non-responding pairs, and the budget spent, donors and batches used (the design ledger).
