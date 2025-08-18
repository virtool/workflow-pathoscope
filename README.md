# workflow-pathoscope

An analysis workflow for detecting known OTUs (viruses) by matching reads against a reference
of known virus genomes.

## Steps

1. **Map default isolates**. Identify potentially present OTUs by mapping against representative (default)
isolates of each OTU.
2. **Build isolate index**. Build a mapping index of all isolates of the previously identified OTUs.
3. **Map all isolates**. Map sample reads against all isolates of OTU candidates.
4. **Map and eliminate subtractions**. Map reads with alignment in previous step against the user-selected subtraction
to eliminate contaminating reads (usually host).
5. **Reassignment**. Use Pathoscope 2.0 to statistically reassign read weight to the most likely reference
genomes of origin. Minimize the impact of multi-mapping and similar reference genomes on the analysis.

## Development

### Container Infrastructure

The project uses Docker containers for development and testing to ensure consistent environments. The development container (`dev` service) automatically starts when running tests and includes all necessary dependencies including Rust toolchain, Python, and bioinformatics tools (bowtie2, HMMER, FastQC, pigz).

The container mounts the project directory and uses persistent volumes for build artifacts (`.venv`, `target`, and UV cache) to speed up subsequent builds. All test commands (`mise run test`, `mise run test:rust`, `mise run test:python`) automatically ensure the development container is running via the `dev:start` dependency. You can also interact with the container directly using `mise run dev` for an interactive shell, or manage it with `mise run dev:stop`, `mise run dev:clean`, and `mise run dev:rebuild`.

### Testing

Run all tests (Rust + Python):

```bash
mise run test
```

Run only Rust unit tests:

```bash
mise run test:rust
```

Run only Python/e2e tests:

```bash
mise run test:python
```

Pass arguments to pytest:

```bash
mise run test:python -- -vv
mise run test:python -- --snapshot-update
mise run test:python -- tests/test_workflow.py::test_map_isolates
```

Run clippy for linting:

```bash
mise run clippy
```

## Contributing

### Commits

All commits must follow the [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0) specification. See the [Virtool development documentation](https://dev.virtool.ca/en/latest/commits_releases.html) for detailed commit guidelines and examples.