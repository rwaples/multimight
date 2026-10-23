# multimight

Turns a [simACE](https://github.com/rwaples/simACE) replicate into an
[EPIMIGHT](https://github.com/BioPsyk/epimight) input pool. Each row is one
proband × trait × kind of relative, written as a gzipped TSV that
`epimight::Pipeline$new(pool = ...)` accepts.

## Run it

The converter needs R with the `arrow` and `data.table` packages:

```bash
Rscript R/simace_to_epimight.R --cohort=rep1/cohort.parquet --pedigree=rep1/pedigree.parquet
```

With [pixi](https://pixi.sh), `pixi run convert` runs the same script in a pinned environment:

```bash
pixi run convert --cohort rep1/cohort.parquet --pedigree rep1/pedigree.parquet --out OUT.tsv.gz
```

`--help` lists every option. The defaults:

| Option | Default | Meaning |
|---|---|---|
| `--out` | `OUT.tsv.gz` | output file |
| `--rels` | `PO,FS,HS,mHS,pHS,Av,1C` | kinds of relatives; `1G` (grandparents) is also available |
| `--traits` | `1,2` | simACE traits to export |
| `--trait-names` | `trait1,trait2` | names written in the `trait` column |
| `--first-birth-year` | `1960` | birth year of the first phenotyped generation |
| `--generation-years` | `30` | years between generations |
| `--death-status` | `2` | `trait_status` for a death before onset: 2 = competing risk, 0 = censored |

The script prints the `relatedness` value to pass to EPIMIGHT's `run_h2` for each kind.

## Output columns

| Column | Meaning |
|---|---|
| `person_id` | simACE id, as text |
| `trait` | trait name |
| `trait_status` | 0 = censored at end of follow-up, 1 = affected, 2 = died before onset |
| `trait_age` | completed years of age at onset, death or censoring |
| `relatives_kind` | the `--rels` code |
| `relatives_n` | relatives of that kind who are phenotyped |
| `relatives_n_trait` | how many of them are affected |
| `birth_year` | from the generation, see `--first-birth-year` |
| `sex` | 0 = female, 1 = male |

## Who counts as a relative

Probands are the phenotyped people: rows of `cohort.parquet` with
`affected1` present. A relative counts only if they are phenotyped too.
Everyone stays in the pool, including people with `relatives_n = 0`,
because EPIMIGHT's population cumulative incidence uses every row.

| Code | Relatives |
|---|---|
| `PO` | mother and father |
| `FS` | full siblings, MZ co-twin left out |
| `HS` | half siblings, maternal and paternal |
| `mHS` | maternal half siblings |
| `pHS` | paternal half siblings |
| `Av` | aunts and uncles: a parent's full sibling or MZ co-twin |
| `1G` | grandparents |
| `1C` | first cousins: an aunt or uncle's children, or a pair sharing two grandparents through half-sibling parents |

`PO`, `Av` and `1G` count only the older relative. The proband is the
child, niece or nephew, or grandchild.

The classification matches
[pedigree-graph](https://github.com/rwaples/pedigree-graph), which simACE
uses. A pair is counted under its closest relationship, so a parent is
never also an aunt or uncle. A parent's MZ co-twin counts as an aunt or
uncle at the nominal relatedness 0.25, although the true value is 0.5.

With two phenotyped generations, no grandparent is phenotyped. That is
why `1G` is not in the default `--rels`. In the first phenotyped
generation, `PO` gives `relatives_n = 0` for the same reason.

## Check it

`check/relative_counts.py` compares every proband's counts with
pedigree-graph's on the same replicate, and exits 1 on any difference:

```bash
pixi run check rep1 OUT.tsv.gz
pixi run check rep1 OUT.tsv.gz early=1 late=2   # with --trait-names early,late
```

## License

MIT
