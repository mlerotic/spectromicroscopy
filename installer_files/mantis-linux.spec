# -*- mode: python -*-

from PyInstaller.utils.hooks import collect_submodules

hidden_imports = collect_submodules('h5py')

block_cipher = None


a = Analysis(['mantis.py'],
             pathex=['/home/watts/dev/mantis/spectromicroscopy'],
             hiddenimports=hidden_imports,
             hookspath=None,
             runtime_hooks=None,
             cipher=block_cipher)

##### include mydir in distribution #######
def extra_datas(mydir):
    def rec_glob(p, files):
        import os
        import glob
        for d in glob.glob(p):
            if os.path.isfile(d):
                files.append(d)
            rec_glob("%s/*" % d, files)
    files = []
    rec_glob("%s/*" % mydir, files)
    extra_datas = []
    for f in files:
        extra_datas.append((f, f, 'DATA'))

    return extra_datas
###########################################

# append the 'data' dir
a.datas += extra_datas('images')
a.datas += extra_datas('file_plugins')
a.datas += extra_datas('TomoCS')

pyz = PYZ(a.pure,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='mantis-linux-3.0.02-x86_64',
          debug=False,
          strip=None,
          upx=True,
          console=False,
          icon='images/Mantis_logo.ico')
