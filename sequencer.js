var fs = require('fs');

if(process.argv.length < 3){
	console.log("Usage:  node sequencer.js  <FASTA filename> <options>");
	console.log("  Application that ingests a FASTA file of DNA segments and attempts");
	console.log("  to reassemble the DNA");
	console.log("Valid options are:");
	console.log("  -v  verbose mode");
	process.exit();
}

var source = process.argv[2];
var fverbose = false;

switch(process.argv[3]){
	case '-v': fverbose = true;
}
const fSep = ':';

/* gDnaSegments is an array of objects that contain the FASTA data from the file
    object: {
	   seq : "string" // the DNA alleles
	   name: "string" // the FASTA segment name
   }
*/
var gDnaSegments = [];

/* gSeedPairs is an array of the initial set of (re)paired DNA sequences from the FASTA file
*/
var gSeedPairs = [];

// readFasta(fname)
//  Function to asynchronously read the data from the FASTA file.  Note that this
//  procedure returns a promise.  Successful read should be handled by the Resolve
//  clause of the promise.
//  call with:
//      fname : the proper string path to the FASTA data file
function readFasta(fname) {
	var p = new Promise(function(resolve,reject){
		fs.readFile(fname, {flag:'r+'}, function(err,data) {
		  if(err) {
			  reject(err);
		  }
		  resolve( processIngest( String(data) ) );
		});
	});
	return p;
}

/* processInject(data)
    This function takes the input stream, and extracts the FASTA data per HUPO-PSI standard
*/
function processIngest(data) {
  var segments = [];
  // Eliminate the Return chars...
  data = data.replace(/\r/g,"");
  
  // if first char is a ';', then remove all chars until first 
  // record start delimiter
  if(data || (data[0] = ';')) {
    var rstart = data.indexOf('>');
	data = (rstart>=0 ? data.substring(rstart, data.length-rstart) : "");
  }
  
  // Split into individual records
  var farray = data.split('>');
  farray.shift();  // eliminate the empty record caused by first element's '>'
  var nRec = farray.length;
  
  // For each FASTA Record, recompose string
  //
  for(var i=0; i<nRec; i++) {
	  var aRec = {};
	  var lines = farray[i].split('\n');
	  aRec.name = lines[0];
	  aRec.seq  = "";
	  for(var s=1; s<lines.length; s++) {
		  aRec.seq += lines[s];
	  }
	  segments.push(aRec);
  }
  return segments;
}


/* overlaps(sumSegment,newSegment)
    This routine performs a string overlap comparison to determine if the
	sequences have a common overlapping set of Alleles.
	Call with:
	   sumSegment : object with running sum of DNA Alleles
	   newSegment : the object with FASTA segment to test for splicing
	   
	Returns:
	   null  : no common overlap encompassing >50% of newSegment
	   -or-
	   {} object with the spliced sum and new segments
*/
function overlaps(sumSegment,newSegment) {
	var a = sumSegment.seq;
	var b = newSegment.seq;
	var lena = a.length;
	var lenb = b.length;
	var half = Math.floor((b.length+1) / 2);
	
	var res = null;
	for(var t=1; t<half; t++){
		var bclip = b.slice(0,lenb-t);
		var m = a.indexOf(bclip);
		
		if(m<0){
			continue;
		}
		if( m===(lena - bclip.length) ) {
			res = {};
			res.name = sumSegment.name + fSep + newSegment.name;
			res.seq  = a + b.slice(lenb-t,lenb);
			break;
		}
	}
	return res;
}


function showCombinations( ds ) {
	var cb = ds.comboPairs;
	console.log("Pair ombinations: " );
	for(var r=0; r<cb.length; r++){
		var s= r+") ";
		for(var c=0; c<cb[r].length; c++){
			var p = cb[r][c];
			s += " [" + p.a + "," + p.b + "]";
		}
		console.log(s);
	}
}


// createSeedPairs()
//  This routine exhaustively walks through all the FASTA record combinations and calculates
// the set of possible initial DNA sequence pairings.
//
function createSeedPairs() {
	var i;
	var row;
	var combos = [];

	// find all the "seedPairs" duple combinations and save
	for(n=0; n<gDnaSegments.length; n++){


	row =[];
		for(i=1; i<gDnaSegments.length; i++){
			var p ={
				a:n,
				b:i
			}
			if(i != n){
				row.push(p);		
			}
		}
		
		for(i=0; i<row.length; i++) {
			var pair = row[i];
		
			var ret = overlaps(gDnaSegments[pair.a], gDnaSegments[pair.b]);
			if(ret != null) {
				gSeedPairs.push(pair);
			}
		}
	}
}


/*
findChain( testCase )
   This is the key routine that recursively tries to append DNA segments onto a running sum.
   If the recursion runs to exhaustion of all segments, testCase returns the Object with the
   sequenced Alleles as well as an ordered string of the FASTA segment-names.
   
   Call with:
      testCase : {
		  dna : {
			  name: "string"  // FASTA segment name(s)
			  seq : "string"  // the DNA sequence
		  },
		  availability : [boolean]  // boolean array indicating if DNA segment is availabe to splice
	  }
	  
	Returns:
	  NULL : if unable to create complete DNA sequence, ordered
	  An object containing the search results
			{
				dna : {},
				availability: []
			}
*/

function findChain( testCase ){
  var moreSegments= function() {
	  var more = false;
	  for(var i=0; i<testCase.availability.length; i++){
		  if(testCase.availability[i]){
			  more = true;
			  break;
		  }
	  }
	  return more;
  };
  
  var next = {};
  
  // No more segments then we're done
  if( !moreSegments() ){
	  return testCase;
  }
  
  // go through all available global segments, see if one will splice
  for(var i=0; i<testCase.availability.length; i++){

	  // Here's an available global DNA segment
	  if(testCase.availability[i]){
		  
		  // try splicing...
		  var next = {
			  dna : overlaps(testCase.dna,gDnaSegments[i]),
			  availability : []
		  }
		  
		  // console.log(testCase.dna.name + " --> tested gDnaSegments["+i +"], result=" + next.dna);
		  
		  if(next.dna != null){
			  // It spliced.  Now mark that global DNA as unavailable
			  for(var p=0; p<testCase.availability.length; p++){
				  next.availability.push(testCase.availability[p]);
			  }
			  next.availability[i] = false;

	  		  // console.log("DNA extended to: " + next.dna.name);
			  
			  // and recurse with new DNA
			  return findChain(next);
		  }
	  }
  }
  // If we get here, then there were no matches
  return null;
}
  
  
function listAlleles(dna){
	
	process.stdout.write("\nRebuilt DNA Sequence from FASTA Segments\n");
	process.stdout.write(  "vv-------------clip here--------------vv\n");
    process.stdout.write(dna.seq);	
	process.stdout.write("\n");
	
}
  
  
function listSeqNames(dna){
	 var fasta = dna.name.split(fSep);
	 var spp= 5;

	process.stdout.write("\nOrdered listing of FASTA segments\n");
	process.stdout.write(  "---------------------------------\n");
	
	for(var i=0; i<fasta.length; i++){
		if( (i>0) && (i%spp==0)){
			process.stdout.write("\n");
		}
		process.stdout.write(fSep + fasta[i]);
	}
	process.stdout.write("\n");
}


console.log("\n****************************************************");
console.log(  "* FASTA Resequencer Application");
console.log(  "* 2017, g. Inn");
console.log( "*");
console.log(  "****************************************************");
console.log("Processing file: " + source);

function analyzeSet(){
	readFasta(source).then( function(segments){
		var theDna = {
			name:"",
			seq:""
		};
		var availability = [];
		var fullChain={};
		
		// Save data from file into global DNA segment array
		gDnaSegments = segments;
		console.log("Read [ " + gDnaSegments.length + " ] DNA sequences from: " + source);


		// Create set of "starting" DNA pairs
		createSeedPairs();
		if(fverbose){
			console.log("There are initially < " + gSeedPairs.length + " > seed pairs of DNA sequences\n");
		}
		
		for(var p=0; p<gSeedPairs.length; p++){
			// Now define the search...
			var testCase = {
				dna : overlaps( gDnaSegments[gSeedPairs[p].a], gDnaSegments[gSeedPairs[p].b]),
				availability : []
			}

			// Create status array of segment availability to splice onto a chain
			for(var i=0; i< gDnaSegments.length; i++){
				testCase.availability.push(true);
			}
			testCase.availability[gSeedPairs[p].a] = false;
			testCase.availability[gSeedPairs[p].b] = false;
			
			if(fverbose){
				process.stdout.write("Building Total DNA from Seed Pair(" + p + "): " + gDnaSegments[gSeedPairs[p].a].name
						+ fSep + gDnaSegments[gSeedPairs[p].b].name );
			}
			fullChain  = findChain(testCase);
			if(fullChain != null){
				if(fverbose){
					process.stdout.write(" ..success!\n");
				}
				break;
			} else {
				if(fverbose){
					process.stdout.write(" fail.\n");
				}
			}
		}	

        if(fullChain != null){
			if(fverbose){
				console.log("\nSequencer was successful to reconstring the DNA chain sequence")
				listSeqNames(fullChain.dna);
			}
			listAlleles(fullChain.dna);
		}
		console.log("\n*** done ***");
		
	},function(error){
		console.log("Error reading Fasta data: " + error);	
	});
	
}


analyzeSet();
