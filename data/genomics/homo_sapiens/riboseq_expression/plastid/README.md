# plastid module test files

The files in this directory are needed to run nf-test on the plastid module/submodules. They were generated using nf-test, too:

- `Homo_sapiens.GRCh38.111_chr20_rois.txt` was generated using the command `nf-core modules test plastid/metagene_generate`.
- `SRX11780887_p_offsets.txt` was generated using the command `nf-core modules test plastid/psite`.

After running the above commands, the files can be extracted from the `work` directory used by nf-test: `.nf-test/tests/<hash>/work/<hash>/<hash>/`
