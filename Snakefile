#!/usr/bin/env python

rule pyrequiremests:
    shell:
        "pip install -r requirements.txt"

rule jlrequirements:
    shell:
        "julia install_packages.jl"


rule get_copra:
    run:
        import wget
        wget.download("http://gregory.org/research/networks/software/copra.jar", out="chaldea")



