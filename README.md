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

## Contributing

### Commits

All commits must follow the [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0) specification.

These standardized commit messages are used to automatically publish releases using [`semantic-release`](https://semantic-release.gitbook.io/semantic-release)
after commits are merged to `main` from successful PRs.

**Example**

```text
feat: add API support for assigning labels to existing samples
```

Descriptive bodies and footers are required where necessary to describe the impact of the commit. Use bullets where appropriate.

Additional Requirements
1. **Write in the imperative**. For example, _"fix bug"_, not _"fixed bug"_ or _"fixes bug"_.
2. **Don't refer to issues or code reviews**. For example, don't write something like this: _"make style changes requested in review"_.
Instead, _"update styles to improve accessibility"_.
3. **Commits are not your personal journal**. For example, don't write something like this: _"got server running again"_
or _"oops. fixed my code smell"_.

From Tim Pope: [A Note About Git Commit Messages](https://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)