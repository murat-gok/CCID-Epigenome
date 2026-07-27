# Upgrade from v0.4.1 to v0.5.0

1. Back up the current project if desired:

```bash
cp -a ~/ccid_epigenome ~/ccid_epigenome_backup_v0.4.1
```

2. Apply the patch:

```bash
cd ~/ccid_epigenome
unzip -o /mnt/c/Users/user/Downloads/ccid_epigenome_patch_v0.5.0.zip -d ~/ccid_epigenome
```

3. Reinstall and test:

```bash
source .venv/bin/activate
python -m pip install -e .
python -m pytest -q
python -c "import ccid; print(ccid.__version__)"
```

Expected version: `0.5.0`; expected tests: `32 passed`.

4. Run the command in `README_O3_VALIDATION_V0.5.0.md`.
