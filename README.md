# 💻 Chem 281 Lab: Version control and GIT

## 🧪 Goal

The goal of this lab is to:

1. Familiarize yourself with **version control workflows**.
2. Learn how to use **git commands**. 
3. Practice using **cheminformatics code** and **different branches**.
4. Fix and merge the **feature branches** into master.

---

## 🗂️ Provided

- A `docker` file to set up the dev environment.
- Cheminformatics code in `scripts/app.py` and associated tests `test/`.

---

## 💻 Setup
```bash
./docker_build.sh # You may need to chmod +x

./docker_run.sh # You may need to chmod +x

python3 script/molecule.py

Total molecules: 17
Molecules passing Lipinski: 0
Molecules with aromatic ring: 0
```

## ✅ Tasks
In the repo you should see there are 3 branches:
```
git branch

bug-fix-1
draw-mols
master
```
master is the branch where we would like all of our changes to be incorporated from the other feature branches. Bug-fix-1 and draw-mols are the feature branches.

### bug-fix-1
Change onto this branch using `git switch bug-fix-1` and try to run `python3 script/molecule.py` and the tests `pytest test/` and you should see that there is a bug! Use `git bisect` to narrow down which commit introduced the bug and record the commit. Then fix the bug by creating a new commit. To wrap up, merge the branch back into master as 1 commit (there are 5+ commits that were added, but we only want 1 to be shown on the master log).

Background about this branch:
It was created to add clustering to the filtered molecules, so we could see which molecules are similar to each other.

### draw-mols
Change onto this branch using `git switch draw-mols` and you should see a bunch of commits. However, we only really care about the first commit that added the molecule drawning functionality. You should use `git cherry-pick` to select this commit and move it onto the master branch.

Background about this branch:
It was created to help visualize the molecules after they have been filtered.

### master
You may have noticed that the tests on master are actually failing, and once again a bad commit is the culprit! Find which commit introduced the failures by using `git bisect` but this this once you have found it, revert the commit using `git revert`.

### Extra time
Now that the molecule drawing code and clustering code is all available on the master branch, can you update the code so that molecules in the same cluster are exported as separate images? This way we can quickly visualize which molecules are in which cluster.
