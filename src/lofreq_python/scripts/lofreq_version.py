#!/usr/bin/env python
"""Dummy script to derive the LoFreq from setup tools
"""

def main():
    """main function
    """

    # see http://stackoverflow.com/questions/2058802/how-can-i-get-the-version-defined-in-setup-py-setuptools-in-my-package
    import pkg_resources  # part of setuptools
    print pkg_resources.require("LoFreq")[0].version
        
if __name__ == "__main__":
    main()
