# -*- coding: utf-8 -*-

import pkgutil, imp, os

# These variables declare the options that each plugin can claim the ability to handle
actions = ['read','write']
data_types = ['spectrum','image','stack','results']

# Go through the directory and try to load each plugin
plugins = []
for m in pkgutil.iter_modules(path=__path__):
    print "Loading file plugin:", m[1], "...",
    try:
        details = imp.find_module(m[1],__path__)
        # could add checks here to enforce presence of required functions in plugin
        plugins.append(imp.load_module(m[1],*details))
        print "("+plugins[-1].title+") Success!"
    except ImportError as e:
        print "prerequisites not satisfied:", e

# Go through set of plugins and assmeble lists of supported file types for each action and data type
supported_filters = dict([a,dict([t,[]] for t in data_types)] for a in actions)
filter_list = dict([a,dict([t,[]] for t in data_types)] for a in actions)
for P in plugins:
    for action in actions:
        for data_type in data_types:
            if data_type in getattr(P,action+'_types'):
                filter_list[action][data_type].append(P.title+' ('+' '.join(P.extension)+')')
                for ext in P.extension:
                    if ext not in supported_filters[action][data_type]:
                        supported_filters[action][data_type].append(ext)
for action in actions:
    for data_type in data_types:
        filter_list[action][data_type] = ['Supported Formats ('+' '.join(supported_filters[action][data_type])+')']+filter_list[action][data_type]
for data_type in data_types:
    filter_list['read'][data_type].append('All files (*.*)')



def load(filename,stack_object=None,plugin=None,selection=None):
    """
    Pass the load command over to the appropriate plugin so that it can import data from the named file.
    """
    if plugin is None:
        plugin = identify(filename)
    if plugin is None:
        return None
    else:
        print "load", filename, "with the", plugin.title, "plugin."
        if stack_object is None:
            return plugin.read(None,filename,selection)
        else:
            plugin.read(stack_object,filename,selection)
            return

def GetFileStructure(filename,plugin=None):
    """
    Use the plugin to skim-read the file and return the structure of the data.
    Returns None if there is only a single data array (i.e. no choices to be made).
    """
    if plugin is None:
        plugin = identify(filename)
    if plugin is None:
        return None
    else:
        print "get info from", filename, "with the", plugin.title, "plugin."
        FileInfo = plugin.GetFileStructure(filename)
        #if FileInfo is not None:
            #print len(FileInfo), len(FileInfo[next(iter(FileInfo))])
            #print FileInfo
        return FileInfo
    
def identify(filename):
    """
    Cycle through plugins until finding one that claims to understand the file format.
    First it tries those claiming corresponding file extensions, followed by all other plugins until an appropriate plugin is found.
    """
    print "Identifying file:", filename, "...",
    ext = os.path.splitext(filename)[1]
    #print "identify", filename, ext
    flag = [True]*len(plugins)
    for i,P in enumerate(plugins):
        if '*'+ext in P.extension:
            if P.identify(filename):
                print "as type:", P.title
                return P
            else:
                flag[i] = False
    for i,P in enumerate(plugins):
        if flag[i]:
            if P.identify(filename):
                print "as type:", P.title
                return P
    print "Error! unknown file type."
    return None

