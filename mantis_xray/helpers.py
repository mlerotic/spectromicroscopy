def resource_path(relative_path):
    import os, sys

    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(__file__))
    return os.path.join(base_path, relative_path)


def check_for_updates(current_version):
    import urllib.request, re
    # Scrape the version string from the PyPI:mantis-xray RSS feed
    try:
        with urllib.request.urlopen("https://pypi.org/rss/project/mantis_xray/releases.xml") as init_file:
            pypi_rss = init_file.read()
        #print(pypi_rss)
        pypi_list = re.findall(r"(?:\<title\>\s*)([\d\.]+)(?:\s*\</title\>)", pypi_rss.decode())[0]
        print("Latest package on PyPI is version {0}".format(pypi_list))
    except:
        pass
    
    # Scrape version string from the code in the github repository
    try:
        with urllib.request.urlopen("https://raw.githubusercontent.com/mlerotic/spectromicroscopy/master/mantis_xray/__init__.py") as init_file:
            github_init = init_file.read()
        #print(github_init)
        github_latest = re.search(r"(?:__version__*\s=*\s)['|\"]+([\d\.]+)", github_init.decode()).group(1)
        print("Current development code is version {0}".format(github_latest))
    except:
        pass
    
    
