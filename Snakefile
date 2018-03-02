#!/usr/bin/env python
rule requiremests:
    shell:
        "pip install -r requirements.txt"
        "julia install_packages.jl"
    run:
        import wget
        copra_url = "http://gregory.org/research/networks/software/copra.jar"
        wget.download(copra_url, out="./chaldea")


