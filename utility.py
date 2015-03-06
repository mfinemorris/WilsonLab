import os, errno, sys, itertools, contextlib

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


    
def alert(Freq = 1500,Dur = 500):
    '''
    On windows os, make beep at Freq Hz for Dur milliseconds.
    default is Freq @ 2500 Hertz and Dur @ 500 ms
    '''
    try:
        from winsound import Beep
        Beep(Freq,Dur)
    except ImportError as e:
        sys.stdout.write('\a')
        sys.stdout.flush()
    except Exception as e:
        pass

   
def get_folders_files(folder):
    dir_content = os.listdir(folder)
    dir_files = []
    dir_folders = []
    for f in dir_content:
        f = os.path.join(folder,f)
        if os.path.isdir(f):
            dir_folders.append(f)
        else:
            dir_files.append(f)
    return dir_folders, dir_files


def display_content(dir_files, data_type=''):
    #get dirname
    print "Folder: ",os.path.dirname(dir_files[0]), '\n'
    for n,i in enumerate(dir_files):
        file_name = os.path.basename(i)
        if data_type in file_name:
            print n,": ",i

def strings_containing(strings, elements):
    
    strings_w_all_elems = []
    for string in strings:
        #print 
        #print string
        num_elems = 0
        for elem in elements:
            #print "\t", elem,
            if elem in string:
                num_elems += 1
           
        if num_elems == len(elements):
            strings_w_all_elems.append(string)
            
    return strings_w_all_elems
def test_strings_containing():
    dummy = ['misc-TBModel-sec132-eL-IP0_9-gnaps1_0.txt.csv',\
    'misc-TBModel-sec1320-eL-IP0_95-gnaps1_0.txt.csv',\
    'misc-TBModel-sec1320-eL-IP0_9-gnaps1_0.txt.csv']

    print strings_containing(dummy, ['0_95','1320'])
    print strings_containing(dummy, ['1320'])
    print strings_containing(dummy, ['0_95'])


def get_files(start, folder_crit=[], file_crit=[], verbose = True):
    folders = []
    curr_dir = os.getcwd()
    folder_num = 0
    for n, (root, dirs, files) in enumerate(os.walk(start)):
        
        #if root does not contain folder criteria, skip it
        if not strings_containing([root], folder_crit):
            continue
        
        if verbose: print "["+str(folder_num)+"]", root
        folder_num += 1
        file_paths = []
        
        goodfiles = strings_containing(files,file_crit)
        for nn,f in enumerate(goodfiles):
            #remove zip files
            if '.zip' in f:
                files.remove(f)
                continue
            file_path = "r'"+os.path.join(curr_dir,root, f).replace('C:','')+"'"
            if verbose: 
                print " ["+str(nn)+"]", file_path
            file_paths.append(file_path.lstrip('r').replace('\'',''))
        if verbose: 
            print
        if len(file_paths) > 0:
            folders.append(file_paths)
    return list(itertools.chain(*folders))


@contextlib.contextmanager
def stdout_redirected(new_stdout):
    '''
    Redirect print statements and other standard output
    to a new standard output new_stdout
    '''
    save_stdout = sys.stdout
    sys.stdout = new_stdout
    try:
        yield None
    finally:
        sys.stdout = save_stdout
