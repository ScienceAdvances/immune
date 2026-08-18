# Contributing

## Development environment

```bash
python -m pip install -e '.[dev,bulk,plot]'
```

Add only the optional scientific backends needed for the change under test.

## Tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

The tests use small in-memory tables and backend fakes for dispatch behavior.
Changes to a mature-library adapter should also receive at least one integration
check against the supported backend version.

## Style

```bash
ruff check src tests examples
```

Public functions should include a concise scientific purpose, parameter
semantics, return type and any count-unit or statistical assumption. New public
functions belong in the appropriate `io`, `pp`, `tl` or `pl` namespace.

## Design principles

1. Keep clone definitions explicit.
2. Keep reads, molecules and cells distinguishable.
3. Preserve original source identifiers.
4. Prefer mature statistical backends over local reimplementation.
5. Return tidy tables or native backend objects.
6. Surface unmatched/ambiguous integration cases.
7. Test the public Scanpy-style API.
