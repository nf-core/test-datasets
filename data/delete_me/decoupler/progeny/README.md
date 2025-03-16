Data obtained with:
``` python
progeny = dc.get_progeny(organism='human', top=1000)
pd.DataFrame(progeny).to_csv("progeny_data.tsv", sep="\t", index=False)
```
