import os
import subprocess


def main():
    source_dir = os.path.abspath(os.path.dirname(__file__))
    # Project root is three levels up
    if os.environ.get("READTHEDOCS_REPOSITORY_PATH"):
        repo_root = os.environ.get("READTHEDOCS_REPOSITORY_PATH")
    else:
        repo_root = os.path.abspath(os.path.join(source_dir, "..", ".."))

    scripts_dir = os.path.relpath(os.path.join(repo_root, "src", "fairmd"), source_dir)
    developer_dir = os.path.relpath(os.path.join(repo_root, "developer"), source_dir)
    auto_dir = "auto_gen"
    template_dir = "_templates/apidoc"
    excluded_patterns = ["*tests*", "*init*"]

    runs = [
        (scripts_dir, "fairmd"),
        (developer_dir, "Developer"),
    ]

    for path, tocname in runs:
        cmd = [
            "sphinx-apidoc",
            "-o",
            auto_dir,
            "--separate",
            "--implicit-namespaces",
            "-d",
            "1",
            "--templatedir",
            template_dir,
            "--tocfile",
            tocname,  # unique per section
            path,
            *excluded_patterns,
        ]
        print("Running command (cwd=docs/src):")
        print(" ".join(cmd))
        subprocess.run(cmd, check=True, cwd=source_dir)
        print("sphinx-apidoc completed successfully for", tocname)


if __name__ == "__main__":
    main()
