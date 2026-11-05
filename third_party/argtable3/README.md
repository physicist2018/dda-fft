# argtable3

Vendored copy of [argtable3](https://github.com/argtable/argtable3), release
**v3.3.1**, taken from the official amalgamation archive
`argtable-v3.3.1-amalgamation.tar.gz` (`dist/argtable3.c`, `dist/argtable3.h`).

- License: BSD-3-Clause, see `LICENSE`.
- Only the two amalgamated files are kept; examples and tests are not.
- Built as C (`cc -std=c11`), linked into the C++ program; the header already
  has `extern "C"` guards.

To update, replace both files with the amalgamation of the new release and
adjust the version above.
