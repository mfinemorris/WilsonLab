import operator
import os
import re
import shutil
from numpy import *
import pandas

debug = True
element_sep = "-"
#period_sub = "o" #may use in future

def preprocess_params(parameters):
    parameters = sorted(parameters.items(), key=operator.itemgetter(1)) #sort dictionary so that items with None will appear first
    params = element_sep.join([key+str(val) if val else key for key,val in parameters]).replace('.','_') #process params
    return params



def make_folderpaths(model, time, parameters={}, path=''):
    params = preprocess_params(parameters)
    time_seconds = str(int(round(time/1000.))) #convert time to seconds and string
    name_parts = [model, "sec"+ time_seconds]
    folder_name = element_sep.join((name_parts+[params] if params else name_parts))
    if path and os.path.isdir(path):
        return os.path.join(path, folder_name)
    else:
        return folder_name



def make_filepaths(info_type, model, time, parameters={}, extension='txt', path=''):
    '''
    Creates a complete path name for a file, where the file name has the following format:

    info_type '-' model '-sec' time '-' parameters '.' extension

    info_type = string indicating what data is being saved. 'voltage' or 'timing' or something else.
    model = the model that was simulated
    time = number of milliseconds that were simulated
    parameters = is a dictionary where any parameter changes can be specified as key:value pairs (key is the name of the parameter!). 
    If the data comes from runs with different values of the same parameter, the value should be None.
    extension = the extension WITHOUT the '.' i.e. 'txt', 'csv', etc.
    path = the path if different from the current working directory
    '''
    #prepare file name components
    params = preprocess_params(parameters)
    time_seconds = str(int(round(time/1000.))) #convert time to seconds and string
    file_name_parts = [info_type, model, "sec"+ time_seconds]

    file_name = element_sep.join((file_name_parts+[params] if params else file_name_parts)) + "." + extension
    if path and os.path.isdir(path):
        return os.path.join(path, file_name)
    else:
        return file_name

def deconstruct_filepaths(file_path, element_seperator=None):
    '''
    Deconstructs a filepath of format:

    info_type '-' model '-sec' time '-' parameters '.' extension

    and returns 

    '''
    
    ## break up the path into its components ##
    #separate the directory from the file name
    base = os.path.basename(file_path)
    directory = os.path.dirname(file_path)
    #separate the file extension and file path
    file_info, delim, ext = base.partition('.')
    #break the file name into its components
    if element_seperator:
        file_info = file_info.split(element_seperator)
    else:
        file_info = file_info.split(element_sep)

    stored_data_type = file_info[0]
    model = file_info[1]
    
    params = {}
    for i in file_info[2:]:
        param_name = re.sub("[^A-z]", "", i).replace("_","") #extract secondary parameter from file name
        param_val = re.sub("[^0-9.]", "", i.replace("_",".")) #extract secondary parameter from file name
        params[param_name] = param_val
    seconds = float(params['sec'])*1000.
    del params['sec']
    
    return stored_data_type, model, seconds, params, ext, directory

def test():
    #test that it works
    print "Folder: ",make_folderpaths('BRSModel', 1000*60*6, parameters={'IP':None, 'gL':2.3, 'gnaps':None}, path = os.getcwd())

    print "Calling make_filepaths('voltage', 'BRSModel', 1000*60*6, parameters={'IP':None, 'gL':2.3, 'gnaps':None}) creates name:"
    print
    path = make_filepaths('voltage', 'BRSModel', 1000*60*6, parameters={'IP':None, 'gL':2.3, 'gnaps':None}, path = os.getcwd())
    print "constructed path: ",path
    print
    path_parts = deconstruct_filepaths(path)
    print path_parts
    print
    print make_filepaths(*path_parts)


def convert_to_new_scheme(file_dir, pattern):

    #get all files in file_dir matching pattern
    dir_contents = os.listdir(file_dir) #get list of all files in file_dir
    matching_files =  [data_file for data_file in dir_contents if re.search(pattern, data_file)] #make list of files matching pattern
    
    for datafile in matching_files:
        #separate the file extension and file path
        name_stem, delim, ext = datafile.partition('.')
        #break the file name into its components
        file_info = name_stem.split("-")
        #format: info_type '-' model '-sec' time '-' parameters '.' extension
        try:
            new_file = "-".join([file_info[0],file_info[2],file_info[3],file_info[1]]+file_info[4:])+delim+ext
        except:
            print "Could not gen new file name for", datafile
            continue
        
        old_path, new_path = os.path.join(file_dir, datafile), os.path.join(file_dir, new_file)
        
        try:
            os.rename(old_path, new_path)
            print "Renamed", old_path, "to", new_path
        except:
            print "Failed to rename", old_path, "to", new_path









if __name__ == "__main__":
    if debug:
        test()
        sys.exit()
    #convert_to_new_scheme(os.path.join(os.getcwd(), 'Results'), 'voltages-IP-TBModel')
    

