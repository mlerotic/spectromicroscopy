# -*- coding: utf-8 -*-

import pkgutil, imp, os

filters = {}
plugins = []
for m in pkgutil.iter_modules(path=__path__):
    print "Loading file plugin:", m[1], "...",
    try:
        details = imp.find_module(m[1],__path__)
        plugins.append(imp.load_module(m[1],*details))
        print "("+plugins[-1].title+") Success!"
    except ImportError as e:
        print "prerequisites not satisfied:", e

def load(filename,plugin=None):
    """
    Use the plugin to read the file and return the data structure passed by the plugin.
    """
    if plugin is None:
        plugin = identify(filename)
    if plugin is None:
        return None
    else:
        print "load", filename, "with the", plugin.title, "plugin."
        return plugin.read(filename)
    
def identify(filename):
    """
    Cycle through plugins until finding one that claims to understand the file format.
    First it tries those claiming corresponding file extensions, followed by all other plugins until an appropriate plugin is found.
    """
    print "Identifying file:", filename
    ext = os.path.splitext(filename)[1]
    #print "identify", filename, ext
    flag = [True]*len(plugins)
    for i,P in enumerate(plugins):
        if '*'+ext in P.extension:
            if P.identify(filename):
                return P
            else:
                flag[i] = False
    for i,P in enumerate(plugins):
        if flag[i]:
            if P.identify(filename):
                return P
    print "Error! unknown file type."
    return None

class File_GUI():
    """
    Ask user to choose file and then use an appropriate plugin to read and return a data structure.
    """
    def __init__(self):
        actions = ['read','write']
        data_types = ['spectrum','image','stack','results']
        self.last_path = dict([a,dict([t,os.getcwd()] for t in data_types)] for a in actions)
        self.last_filter = dict([a,dict([t,0] for t in data_types)] for a in actions)
        self.supported_filters = dict([a,dict([t,[]] for t in data_types)] for a in actions)
        self.filter_list = dict([a,dict([t,[]] for t in data_types)] for a in actions)
        for P in plugins:
            for action in actions:
                for data_type in data_types:
                    if data_type in getattr(P,action+'_types'):
                        self.filter_list[action][data_type].append(P.title+' ('+' '.join(P.extension)+')')
                        for ext in P.extension:
                            if ext not in self.supported_filters[action][data_type]:
                                self.supported_filters[action][data_type].append(ext)
        for action in actions:
            for data_type in data_types:
                self.filter_list[action][data_type] = ['Supported Formats ('+' '.join(self.supported_filters[action][data_type])+')']+self.filter_list[action][data_type]
        for data_type in data_types:
            self.filter_list['read'][data_type].append('All files (*.*)')
    
    def SelectFile(self,action,data_type):
        from PyQt4 import QtGui
        dlg=QtGui.QFileDialog(None)
        dlg.setWindowTitle('Choose File')
        dlg.setViewMode(QtGui.QFileDialog.Detail)
        dlg.setDirectory(self.last_path[action][data_type])
        dlg.setFilters(self.filter_list[action][data_type])
        dlg.selectFilter(self.filter_list[action][data_type][self.last_filter[action][data_type]])
        if dlg.exec_(): #if not cancelled
            self.last_path[action][data_type] = os.path.split(str(dlg.selectedFiles()[0]))[0]
            chosen_plugin = None
            for i,filt in enumerate(self.filter_list[action][data_type][1:-1]):
                if filt==dlg.selectedFilter():
                    chosen_plugin = plugins[i]
                    break
            if chosen_plugin is not None:
                self.last_filter[action][data_type] = i+1
            return (str(dlg.selectedFiles()[0]),chosen_plugin)
        else:
            print "cancelled"
            return (None,None)

