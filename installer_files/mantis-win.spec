# -*- mode: python -*-

block_cipher = None


a = Analysis(['mantis_qt.py'],
             pathex=['C:\\projects\\Mantis'],
             hiddenimports=['scipy.special._ufuncs_cxx','h5py.h5ac'],
             hookspath=None,
             runtime_hooks=None)

for d in a.datas:
    if 'pyconfig' in d[0]: 
        a.datas.remove(d)
        break

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

pyz = PYZ(a.pure)

exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='mantis.exe',
          debug=False,
          strip=None,
          upx=True,
          console=False,
          icon='images/Mantis_logo.ico' )
