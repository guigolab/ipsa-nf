import nextflow.Nextflow

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
		  def path = Nextflow.file(list[2])
		  def type = list[3]
		  def view = list[4]
		  def readType = list[5]
		  def readStrand = list[6]
		  def readLength = list[7]
		  bams << [ mergeId, id, path, type, view, readType, readStrand, readLength ]
		}
		return bams
	}
}