#!/usr/bin/env groovy
// Exercise the actual labelled selection section, without loading a helper.
String source = new File('main.nf').getText('UTF-8')
String begin = '// BEGIN INPUT SELECTION'
String end = '// END INPUT SELECTION'
assert source.count(begin) == 1
assert source.count(end) == 1
String selection = source.substring(source.indexOf(begin) + begin.length(), source.indexOf(end))
Closure resolve = { Map options, String table ->
    Binding binding = new Binding()
    binding.setVariable('params', options)
    binding.setVariable('file', { Object... arguments ->
        Object requested = arguments.find { !(it instanceof Map) }
        assert requested.toString() == options.sequence_inputs
        if (table == null) throw new IllegalArgumentException('Cannot read --sequence_inputs.')
        [text: table]
    })
    binding.setVariable('error', { Object message -> throw new IllegalArgumentException(message.toString()) })
    new GroovyShell(binding).evaluate(selection, 'main-input-selection.groovy')
    [inputMode: binding.getVariable('input_mode'), treeSource: binding.getVariable('tree_source')]
}
String focal = 'species_id\tinput_file\tinput_type\trole\nA\ta.fa\tgenome\tfocal\nB\tb.fa\tgenome\tfocal\n'
String genomes = focal + 'O\to.fa\tgenome\toutgroup\n'
String proteomes = genomes.replace('genome', 'proteome')
String mixed = focal + 'O\to.faa\tproteome\toutgroup\n'
int cases = 0
Closure check = { Map options, String table, String mode, String tree ->
    Map result = resolve(options, table)
    assert result.inputMode == mode
    assert result.treeSource == tree
    cases++
}
Closure reject = { Map options, String table, String message ->
    boolean rejected = false
    try { resolve(options, table) }
    catch (IllegalArgumentException e) {
        assert e.message.contains(message): e.message
        rejected = true
    }
    assert rejected: "Expected rejection: ${options}"
    cases++
}
check([sequence_inputs:'s.tsv', annotation:'a.gff'], genomes, 'annotation', 'inferred')
check([sequence_inputs:'s.tsv', annotation:'a.gff', species_tree:'t.nwk'], focal, 'annotation', 'supplied')
check([protein_fasta:'p.faa', id_species:'m.tsv', sequence_inputs:'s.tsv'], genomes, 'protein', 'inferred')
check([protein_fasta:'p.faa', id_species:'m.tsv', species_tree:'t.nwk'], null, 'protein', 'supplied')
check([protein_fasta:'p.faa', id_species:'m.tsv', sequence_inputs:'s.tsv'], proteomes, 'protein', 'inferred')
check([sequence_inputs:'s.tsv', annotation:'a.gff'], mixed, 'annotation', 'inferred')
check([protein_fasta:'p.faa', id_species:'m.tsv', sequence_inputs:'s.tsv'], mixed, 'protein', 'inferred')
check([sequence_inputs:'s.tsv', annotation:'a.gff', species_tree:'t.nwk'], genomes, 'annotation', 'supplied')
check([protein_fasta:'p.faa', id_species:'m.tsv', sequence_inputs:'s.tsv', species_tree:'t.nwk'], genomes, 'protein', 'supplied')
check([protein_fasta:'p.faa', id_species:'m.tsv', sequence_inputs:'unused.tsv', species_tree:'t.nwk'], null, 'protein', 'supplied')
check([sequence_inputs:'s.tsv', annotation:'a.gff', species_tree:'t.nwk'], focal + 'O\t\tinvalid\toutgroup\n', 'annotation', 'supplied')
check([sequence_inputs:'s.tsv', annotation:'a.gff', species_tree:'t.nwk'], genomes + 'O2\tx.fa\tgenome\toutgroup\n', 'annotation', 'supplied')
check([sequence_inputs:'s.tsv', annotation:'a.gff', protein_fasta:'incomplete.faa'], genomes, 'annotation', 'inferred')
check([protein_fasta:'p.faa', id_species:'m.tsv', annotation:'incomplete.gff', species_tree:'t.nwk'], null, 'protein', 'supplied')
check([protein_fasta:'p.faa', id_species:'m.tsv', annotation:'a.gff', sequence_inputs:'s.tsv'], proteomes, 'protein', 'inferred')
reject([:], null, 'No complete analysis')
reject([protein_fasta:'p.faa', species_tree:'t.nwk'], null, 'No complete analysis')
reject([sequence_inputs:'s.tsv', annotation:'a.gff', protein_fasta:'p.faa', id_species:'m.tsv'], genomes, 'Both analysis')
reject([sequence_inputs:'s.tsv', annotation:'a.gff'], focal, 'Without --species_tree')
reject([protein_fasta:'p.faa', id_species:'m.tsv'], null, 'Without --species_tree')
reject([input_mode:'protein'], null, '--input_mode has been removed')
reject([species_inputs:'old.tsv'], null, '--species_inputs has been removed')
reject([sequence_inputs:'s.tsv'], 'bad\theader\n', 'requires unique columns')
reject([sequence_inputs:'s.tsv'], focal + 'A\tc.fa\tgenome\tfocal\n', 'Duplicate species_id')
reject([sequence_inputs:'s.tsv'], genomes.replace('genome', 'unknown'), 'Unsupported input_type')
reject([sequence_inputs:'s.tsv'], focal + 'O\to.fa\tgenome\tother\n', 'Unsupported sequence-input role')
reject([sequence_inputs:'s.tsv', annotation:'a.gff'], genomes + 'O2\to2.fa\tgenome\toutgroup\n', 'exactly one auxiliary outgroup')
reject([sequence_inputs:'s.tsv', annotation:'a.gff'], proteomes, 'No complete analysis')
println "Passed ${cases} inline input-selection tests."
