import nextflow.Nextflow
import java.nio.file.Path

/**
 *  Helper class to read input files and metadata from a TAB separated file
 *
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
 */
class IPSA {

	/**
	* Define input fields position:
	* 1. id
	* 2. path
	* 3. read type
	* 4. read strand
	*/
	class Fields {
		static final int ID          = 0
		static final int PATH        = 1
		static final int READ_TYPE   = 2
		static final int READ_STRAND = 3
	}

	/**
	* Parse the index TSV file as defined in the Fields class:
	*/
	static parseIndexFile(def f) {
		def bams = []
		f.readLines().each { line ->
			def list = line.split()
		  def id = list[Fields.ID]
		  def path = resolveFile(list[Fields.PATH], f)
		  def readType = list[Fields.READ_TYPE]
		  def readStrand = list[Fields.READ_STRAND]
		  bams << [ id, path, readType, readStrand ]
		}
		return bams
	}

	/*
	 * Given a string path resolve it against the index file location.
	 * Params: 
	 * - str: a string value represting the file pah to be resolved
	 * - index: path location against which relative paths need to be resolved 
	 */
	static resolveFile( String str, Path index ) {
	  if( str.startsWith('/') || str =~ /^[\w\d]*:\// ) {
	    return Nextflow.file(str)
	  }
	  else if( index instanceof Path ) {
	    return index.parent.resolve(str)
	  }
	  else {
	    return Nextflow.file(str) 
	  }
	} 
}