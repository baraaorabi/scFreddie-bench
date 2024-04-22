git clone https://github.com/baraaorabi/scFreddie-bench
cd scFreddie-bench
git submodule update --init --recursive
mamba env create -f env.yaml
mamba activate sc_bench
snakemake -j96 --use-conda --conda-create-envs-only