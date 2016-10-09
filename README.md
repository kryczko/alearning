# alearning (Aluminum learning)

a learning thus far is a collection of scripts with the goal of machine learning DFT calculations of aluminum based systems. This includes aluminum doped and aluminum systems with defects.

Thus far the repo contains 2 scripts:

genfccal.py -> generates an initial perfect crystal of Al

gen_train_data.py -> generates configurations of doped-Al and defect-Al crystals depending on user configuration.

Example:
```bash
$ python gen_train_data.py -de 1 -do Ti 5 O 3 -n 100
```

This generates 100 configurations of Al crystals with one defect, 5 Ti and 3 O dopant atoms.
