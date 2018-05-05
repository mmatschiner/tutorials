## this script takes a batch of (simple) nexus files and concatenates them into a single file of nexus, phylip, or fasta format.

# Define the valid formats.
valid_formats = ["nexus","phylip","fasta"]
format_raise_string = "A file format ("
valid_formats.each {|f| format_raise_string << "#{f}, "}
format_raise_string.chomp!(", ")
format_raise_string << ") must be specified with the option \"-f\"!"

# Define the help string.
help_string = ""
help_string << "\n"
help_string << "concatenate.rb\n"
help_string << "\n"
help_string << "Available options:\n"
help_string << "  Option   Value                    Comment\n"
help_string << "  -i       input file names       | Full or relative path plus file name (wildcards accepted)\n"
help_string << "  -o       output file name       | Full or relative path plus file name\n"
help_string << "  -f       format                 | Output format, can be "
valid_formats[0..-2].each {|f| help_string << "#{f}, "}
help_string << "or #{valid_formats[-1]}\n"
help_string << "\n"

# Read the arguments.
if ARGV == [] or ["-h","--help","-help"].include?(ARGV[0].downcase)
    puts help_string
    exit
end

# Read the specified format and remove the -f option from the ARGV array.
if ARGV.include?("-f")
    raise format_raise_string if ARGV[ARGV.index("-f")+1] == nil
    format = ARGV[ARGV.index("-f")+1].downcase
    unless valid_formats.include?(format)
        puts format_raise_string
        exit 1
    end
else
    puts format_raise_string
    exit
end
ARGV[ARGV.index("-f")+1] = nil
ARGV[ARGV.index("-f")] = nil
ARGV.compact!

# Read the specified output file name.
if ARGV.include?("-o")
    if ARGV[ARGV.index("-o")+1] == nil
        puts "ERROR: Please specify an output file name with the option \"-o\"!"
        exit
    end
    output_file_name = ARGV[ARGV.index("-o")+1]
    ARGV[ARGV.index("-o")+1] = nil
    ARGV[ARGV.index("-o")] = nil
    ARGV.compact!
else
    output_file_name == nil
end

# Read the specified input file names.
if ARGV.include?("-i")
    input_file_names = ARGV[ARGV.index("-i")+1..-1]
else
    puts "ERROR: At least one input file name must be given with option \"-i\"!"
    exit 1
end

# Initiate arrays for output ids and output sequences.
ids_per_file = []
seqs_per_file = []

input_file_names.each do |a|
    print "INFO: Reading #{a}..." unless output_file_name == nil
    file_ids = []
    file_seqs = []
    file = File.open(a)
    lines = file.readlines
    unless lines[0].strip.downcase == "#nexus"
        puts "ERROR: File is not in NEXUS format: #{a}!"
        exit 1
    end
    in_matrix = false
    lines.each do |l|
        l.strip!
        if l.downcase == "matrix"
            in_matrix = true
        elsif l == ";"
            in_matrix = false
        elsif l != "" and in_matrix
            file_ids << l.split[0]
            file_seqs << l.split[1]
        end
    end
    file_seqs.each do |s|
        if s.size != file_seqs[0].size
            puts "ERROR! Two sequences have different lengths!"
            exit 1
        end
    end
    if file_ids.uniq != file_ids
        puts "ERROR! Not all IDs are unique in file #{a}!"
        exit 1
    end
    ids_per_file << file_ids
    seqs_per_file << file_seqs
    puts " done." unless output_file_name == nil
end

# Check if all id arrays are identical.
all_id_arrays_identical = true
ids_per_file.each do |ids|
    if ids != ids_per_file[0]
        all_id_arrays_identical = false
    end
end

# Get an array of only the unique ids.
if all_id_arrays_identical
    unique_ids = ids_per_file[0].dup
else
    unique_ids = []
    ids_per_file.each do |ids|
        ids.each do |id|
            unique_ids << id unless unique_ids.include?(id)
        end
    end
    unique_ids.sort!
end

# Prepare an array of concatenated sequences.
concatenated_seqs = []
unique_ids.each do |uid|
    concatenated_seq = ""
    ids_per_file.size.times do |x|
        if ids_per_file[x].include?(uid)
            concatenated_seq << seqs_per_file[x][ids_per_file[x].index(uid)]
        else
            seqs_per_file[x][0].size.times do
                concatenated_seq << "-"
            end
        end
    end
    concatenated_seqs << concatenated_seq
end

# Determine the maximum id length.
max_id_length = 0
unique_ids.each do |uid|
    max_id_length = uid.size if uid.size > max_id_length
end

out_string = ""
if format == "phylip"
    out_string << "#{unique_ids.size} #{concatenated_seqs[0].size}\n"
    unique_ids.size.times do |x|
        out_string << "#{unique_ids[x].ljust(max_id_length)} #{concatenated_seqs[x]}\n"
    end
elsif format == "nexus"
    out_string << "#nexus\n"
    out_string << "\n"
    out_string << "begin data;\n"
    out_string << "  dimensions ntax=#{unique_ids.size} nchar=#{concatenated_seqs[0].size};\n"
    out_string << "  format datatype=dna gap=- missing=?;\n"
    out_string << "  matrix\n"
    unique_ids.size.times do |x|
        out_string << "  #{unique_ids[x].ljust(max_id_length)} #{concatenated_seqs[x]}\n"
    end
    out_string << "  ;\n"
    out_string << "end;\n"
elsif format == "fasta"
    unique_ids.size.times do |x|
        out_string << ">#{unique_ids[x]}\n"
        seq_to_write = concatenated_seqs[x]
        while seq_to_write.size > 0 do
            out_string << "#{seq_to_write.slice!(0..59)}\n"
        end
    end
else
    puts "ERROR: Unrecognized format #{format}!"
    exit 1
end
if output_file_name == nil
    puts out_string
else
    File.new(output_file_name,"w").puts out_string
    puts "INFO: Wrote file #{output_file_name}."
end

concatenated_seqs.size.times do |x|
    puts "WARNING: No sequence information available for #{unique_ids[x]}!" if concatenated_seqs[x].gsub("-","").length == 0
end
