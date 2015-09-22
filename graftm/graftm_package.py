import json
import os
import shutil

class InsufficientGraftMPackageException(Exception): pass

class GraftMPackage:
    '''Package to represent a GraftM package. To start using a package, run

    pkg = GraftMPackage.acquire('/path/to/a.gpkg')

    and then retrieve values with e.g.

    pkg.alignment_hmm_path()
      #=> '/path/to/a.gpkg/the.hmm'
    '''

    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    # The key names are unlikely to change across package format versions,
    # so store them here in the superclass
    DIAMOND_DATABASE_KEY = 'diamond_database'
    VERSION_KEY = 'graftm_package_version'
    ALIGNMENT_HMM_KEY = 'align_hmm'
    SEARCH_HMM_KEY = "search_hmms"
    REFERENCE_PACKAGE_KEY = "refpkg"
    HMM_TRUSTED_CUTOFF_KEY = "trusted_cutoff"
    RANGE_KEY = "range"
    _CONTENTS_FILE_NAME = 'CONTENTS.json'

    _CURRENT_VERSION = 2

    _REQUIRED_KEYS = {'2': [
                     DIAMOND_DATABASE_KEY,
                     VERSION_KEY,
                     ALIGNMENT_HMM_KEY,
                     SEARCH_HMM_KEY,
                     REFERENCE_PACKAGE_KEY,
                     HMM_TRUSTED_CUTOFF_KEY,
                     RANGE_KEY
                     ]}


    @staticmethod
    def acquire(graftm_package_path):
        '''Acquire a new graftm Package

        Parameters
        ----------
        graftm_output_path: str
            path to base directory of graftm
        '''
        pkg = GraftMPackageVersion2()


        pkg._base_directory = graftm_package_path
        pkg._contents_hash = json.load(
                                       open(
                                            os.path.join(
                                                         graftm_package_path,
                                                         GraftMPackage._CONTENTS_FILE_NAME
                                                         ),
                                             )
                                       )

        # check we are at current version otherwise choke
        pkg.check_universal_keys(2)
        pkg.check_required_keys(GraftMPackageVersion2._REQUIRED_KEYS)
        return pkg

    def check_universal_keys(self, version):
        h = self._contents_hash
        try:
            v = h[GraftMPackage.VERSION_KEY]
        except KeyError:
            raise InsufficientGraftMPackageException("No version information in graftm package")
        if v != version:
            import IPython ; IPython.embed()
            raise InsufficientGraftMPackageException("Bad version: %s" % v)

    def check_required_keys(self, required_keys):
        '''raise InsufficientGraftMPackageException if this package does not
        conform to the standard of the given package'''
        h = self._contents_hash
        for key in required_keys:
            if key not in h:
                raise InsufficientGraftMPackageException("package missing key %s" % key)

    def __getitem__(self, key):
        '''Return the value of the given key from the contents file'''
        return self.contents_hash[key]

class GraftMPackageVersion2(GraftMPackage):
    version = 2

    _REQUIRED_KEYS = [
                     #GraftMPackage.DIAMOND_DATABASE_KEY, #not required for nucleotide packages
                     GraftMPackage.VERSION_KEY,
                     GraftMPackage.ALIGNMENT_HMM_KEY,
                     GraftMPackage.SEARCH_HMM_KEY,
                     GraftMPackage.REFERENCE_PACKAGE_KEY,
                     GraftMPackage.HMM_TRUSTED_CUTOFF_KEY
                     ]

    def diamond_database_path(self):
        if self._contents_hash[GraftMPackage.DIAMOND_DATABASE_KEY]:
            return os.path.join(self._base_directory,
                                self._contents_hash[GraftMPackage.DIAMOND_DATABASE_KEY])
        else:
            return None


    def search_hmm_paths(self):
        return [os.path.join(self._base_directory, x) for x in
                self._contents_hash[GraftMPackage.SEARCH_HMM_KEY]]

    def alignment_hmm_path(self):
        return os.path.join(self._base_directory,
                            self._contents_hash[GraftMPackage.ALIGNMENT_HMM_KEY])

    def reference_package_path(self):
        return os.path.join(self._base_directory,
                            self._contents_hash[GraftMPackage.REFERENCE_PACKAGE_KEY])

    def use_hmm_trusted_cutoff(self):
        return self._contents_hash[GraftMPackage.HMM_TRUSTED_CUTOFF_KEY]

    def maximum_range(self):
        return self._contents_hash[GraftMPackage.RANGE_KEY]
    
    def _refpkg_contents(self):
        return json.loads(open(os.path.join(self.reference_package_path(), 'CONTENTS.json')).read())
    
    def taxtastic_seqinfo_path(self):
        return os.path.join(self.reference_package_path(),
                            self._refpkg_contents()['files']['seq_info'])
        
    def taxtastic_taxonomy_path(self):
        return os.path.join(self.reference_package_path(),
                            self._refpkg_contents()['files']['taxonomy'])
        
    @staticmethod
    def compile(output_package_path, refpkg_path, hmm_path, diamond_database_file, max_range, trusted_cutoff=False):
        '''Create a new GraftM package with the given inputs. Any files
        specified as parameters are copied into the final package so can
        be removed after calling this function.
        
        Parameters
        ----------
        output_package_path: str
            path to the package being created (must not exist)
        refpkg_path: str
            path to pplacer reference package
        hmm_path: str
            path to the search and align HMM
        diamond_database_file: str
            path to diamond DB file, or None for nucleotide packages
        max_rage: str
            as per maximum_range()
        trusted_cutoff: boolean
            set TC in search HMM
            
        Returns
        -------
        Nothing
        '''
        
        if os.path.exists(output_package_path): 
            raise Exception("Not writing new GraftM package to already existing file/directory with name %s" % output_package_path)
        os.mkdir(output_package_path)
        
        hmm_file_in_gpkg = os.path.basename(hmm_path)
        shutil.copyfile(hmm_path, os.path.join(output_package_path, hmm_file_in_gpkg))
        
        if diamond_database_file:
            diamond_database_file_in_gpkg = os.path.basename(diamond_database_file)
            shutil.copyfile(diamond_database_file, os.path.join(output_package_path, diamond_database_file_in_gpkg))
        
        refpkg_in_gpkg = os.path.basename(refpkg_path)
        shutil.copytree(refpkg_path, os.path.join(output_package_path, refpkg_in_gpkg))
        
        contents = {GraftMPackage.VERSION_KEY: GraftMPackageVersion2.version,
                    GraftMPackage.ALIGNMENT_HMM_KEY: hmm_file_in_gpkg,
                    GraftMPackage.SEARCH_HMM_KEY: [hmm_file_in_gpkg],
                    GraftMPackage.REFERENCE_PACKAGE_KEY: refpkg_in_gpkg,
                    GraftMPackage.HMM_TRUSTED_CUTOFF_KEY: trusted_cutoff,
                    GraftMPackage.RANGE_KEY: max_range}
        if diamond_database_file:
            contents[GraftMPackage.DIAMOND_DATABASE_KEY] = diamond_database_file_in_gpkg
        
        json.dump(contents, open(os.path.join(output_package_path, GraftMPackage._CONTENTS_FILE_NAME), 'w'))
    
    

