import nextflow.Nextflow
import java.nio.file.Path

/**
 *  Helper class to read input files and metadata from a TAB separated file
 *
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
 */
class TsvIndexFile {
	
	/**
	*	Parse a TSV file with the following format:
	* 1. merge id
	* 2. run id
	* 3. input file path
	* 4. file type
	* 5. view
	* 6. read type
	* 7. read strand
	* 7. read length
	*/
	static parse(def f) {
		def bams = []
		f.readLines().each { line ->
			def list = line.split()
		  def mergeId = list[0]
		  def id = list[1]
		  def path = resolveFile(list[2], f)
		  def type = list[3]
		  def view = list[4]
		  def readType = list[5]
		  def readStrand = list[6]
		  bams << [ mergeId, id, path, type, view, readType, readStrand ]
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