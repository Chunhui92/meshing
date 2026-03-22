# AGENTS.md

## Purpose

This repository is for meshing, geometric reconstruction, and numerical algorithm work.
Default expectation: research first, design second, implementation third.

Do not jump straight into production code for nontrivial geometry/MIR tasks unless the user explicitly asks for a quick prototype.

## Working Rules

For algorithm-heavy tasks:

1. Clarify the problem, constraints, and success criteria.
2. Review relevant papers, classic methods, and practical implementations when the choice of algorithm is still open.
3. Compare candidate approaches before coding.
4. Define the code structure before editing.
5. Implement incrementally.
6. Run code and validate on representative data.
7. Summarize observed results, limitations, and next steps.

For this repo, do not treat pseudocode or an unrun implementation as "done".

## Codebase Structure

Keep the repository organized like this:

- root:
  - `main.py` only as the main CLI entrypoint
  - project docs such as `MIR_IMPLEMENTATION.md`
- `src/`:
  - core implementation modules
  - geometry utilities
  - reconstruction logic
  - mesh/IO readers
- `test/`:
  - validation scripts
  - VTK analysis scripts
  - visualization helpers

Do not leave new core logic in the repository root.

## MIR-Specific Guidance

The current MIR implementation is a first-stage prototype, not a full production MIR system.

What is implemented:

- lightweight Silo mixed-material reading
- node resampling from cell fractions
- shared-point caching
- `zoo / plane / fast` reconstruction modes
- VTK export
- quantitative VTK analysis and visualization scripts in `test/`

What is not implemented yet:

- full `>2` material zoo reconstruction
- complete VisIt-style zoo case tables
- Isovolume fallback
- iterative volume-fraction correction
- topology-clean final watertight multi-material output

So when extending MIR code:

- distinguish clearly between the target architecture and the current implementation
- do not describe the current code as a full Zoo MIR implementation
- preserve deterministic shared-point behavior
- treat non-manifold interfaces, degeneracies, and volume error as first-class issues

## Implementation Expectations

When editing code:

- preserve numerical robustness over brevity
- keep geometry kernels separate from CLI and test code
- handle degeneracies explicitly
- keep behavior deterministic when possible
- avoid unnecessary dependencies
- prefer small, focused functions and modules

If a method is fragile, either add safeguards or document the limitation in the same round.

## Validation Requirements

Validation is mandatory for geometry/meshing work.

At minimum, check:

- the code runs
- outputs can be opened correctly
- quantitative metrics are reported
- topology-related checks are performed when relevant

For MIR work, use the scripts in `test/`:

- `test/analyze_vtk.py`
  - verifies that generated VTK files are readable
  - reports cell types, material volumes, boundary edges, and non-manifold edges
- `test/visualize_vtk.py`
  - reads VTK output and visualizes full meshes or per-material interfaces

Do not claim success from compilation alone.

## Documentation

If code structure, algorithm behavior, execution commands, or validation results change, update the relevant Markdown files in the same round.

For MIR changes, `MIR_IMPLEMENTATION.md` should stay aligned with:

- current code structure
- what is actually implemented
- measured validation results
- known limitations

## Anti-Patterns

Avoid these unless explicitly requested:

- coding before comparing algorithm choices
- presenting research ideas as if they are implemented
- ignoring tolerances, degeneracies, or topology checks
- leaving validation scripts or docs stale after code changes
