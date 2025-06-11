rule save_config:
    input:
        config_file="config/config.yml"
    output:
        saved_config=f"results/{run}/config.yml"
    shell:
        "cp {input.config_file} {output.saved_config}"
