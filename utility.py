import os, errno

def mkdir_p(path):
    '''
    This function creates a folder at the end of the specified path, unless the folder already exsists. 
    '''
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass #do nothing if the error occurs because the path already exists
        else: raise #re-raises the error


    
