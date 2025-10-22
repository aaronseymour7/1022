#!/usr/bin/env python3

from hfrpkg.multi.grab_unique import grab_unique_coms
from hfrpkg.multi.run_unique import run_jobs

def main():
    grab_unique_coms()

    run_jobs()
def main_cli():
    main()

if __name__ == "__main__":
    main_cli()
