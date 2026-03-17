---
name: Bug report
about: Report something that isn't working as expected
title: "[BUG] "
labels: bug
assignees: ''

---

## Description
Briefly describe what isn't working.

## Expected behavior
What should happen instead?

## Actual behavior
What actually happens? Include error messages if any.

## Steps to reproduce
How can someone reproduce this bug? Be as specific as possible:

1. Load data (example/custom)
2. Navigate to...
3. Click on...
4. Observe error

## Environment
Please provide:
- **R version:** (run `R.version$version.string`)
- **OS:** Windows / macOS / Linux
- **Browser:** (if dashboard accessed remotely)
- **Package versions:** 
  - DESeq2: (run `packageVersion("DESeq2")`)
  - Shiny: (run `packageVersion("shiny")`)
  - ggplot2: (run `packageVersion("ggplot2")`)

## Screenshots or error output
If applicable, paste console errors, warnings, or screenshots here:

```
Paste error message or stack trace here
```

## Data
- Are you using the included example dataset (GSE157234 / Shemer et al. 2020)?
- If custom: what is the approximate size and structure?
- Can the bug be reproduced with the example data?

## Additional context
Anything else that might help debug this issue?

---
**Thank you for reporting!** 🙏 This helps make the project more robust.
