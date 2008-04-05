def outoptionsnobck()
  return "-tis -suf -bwt -lcp -des"
end

def outoptions()
  return outoptionsnobck() + " -bck"
end

def trials()
  return "-scantrials 10 -multicharcmptrials 1000"
end

def checksfx(parts,pl,withsmap,sat,cmp,filelist)
  extra=withsmap
  if cmp
    extra=extra + " -cmpcharbychar"
  end
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run_test "#{$bin}gt suffixerator -v -parts #{parts} -pl #{pl} " +
           "#{extra} #{outoptions()} -indexname sfx -db " + filearg
  run_test "#{$bin}gt dev sfxmap #{trials()} #{outoptions()} -v sfx",
           :maxtime => 600
end

def flattenfilelist(filelist)
  s=""
  filelist.each do |f|
    s += "#{$testdata}#{f} "
  end
  return s
end

def checkbwt(filelist)
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run_test "#{$bin}gt suffixerator -pl #{outoptions()} -indexname sfx -db " +
           flattenfilelist(filelist)
end

def runsfxfail(args)
  Name "gt suffixerator failure"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis " + args,:retval => 1
  end
end

allfiles = ["RandomN.fna","Random.fna","Atinsert.fna",
            "TTT-small.fna","trna_glutamine.fna",
            "Atinsert.fna","Random-Small.fna","Duplicate.fna"]

alldir = ["fwd","cpl","rev","rcl"]

# put the tests with paircmp, maxpair, patternmatch, into a file gt_idxmatch

Name "gt suffixerator paircmp"
Keywords "gt_suffixerator"
Test do
  run_test "#{$bin}gt dev paircmp -a ac 11"
end

Name "gt suffixerator maxpairs"
Keywords "gt_suffixerator"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           "-indexname sfx -dna -tis -suf -lcp -pl"
  run_test "#{$bin}gt dev maxpairs -l 8 -ii sfx"
  run "grep -v '^#' #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}maxpairs-8-Atinsert.txt"
  run_test "#{$bin}gt dev maxpairs -scan -l 8 -ii sfx"
  run "grep -v '^#' #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}maxpairs-8-Atinsert.txt"
  run_test "#{$bin}gt dev maxpairs -samples 40 -l 6 -ii sfx",:maxtime => 600
end

Name "gt suffixerator patternmatch"
Keywords "gt_suffixerator"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           "-indexname sfx -dna -bck -suf -tis -pl"
  run_test "#{$bin}gt dev patternmatch -samples 10000 -minpl 10 -maxpl 15 " +
           " -bck -imm -ii sfx"
end

alldir.each do |dir|
  Name "gt suffixerator #{dir}"
  Keywords "gt_suffixerator"
  Test do
     run_test "#{$bin}gt suffixerator -dir #{dir} -tis -suf -bwt -lcp " +
              "-indexname sfx -pl -db " + 
         flattenfilelist(allfiles)
     run_test "#{$bin}gt suffixerator -storespecialcodes -dir #{dir} -tis -suf -lcp " +
              "-indexname sfx -pl -db " +
         flattenfilelist(allfiles)
     run_test "#{$bin}gt suffixerator -tis -bwt -lcp -pl -ii sfx"
  end
end

runsfxfail "-indexname sfx -db /nothing"
runsfxfail "-indexname /nothing/sfx -db #{$testdata}TTT-small.fna"
runsfxfail "-smap /nothing -db #{$testdata}TTT-small.fna"
runsfxfail "-dna -db #{$testdata}sw100K1.fsa"
runsfxfail "-protein -dir cpl -db #{$testdata}sw100K1.fsa"
runsfxfail "-dna -db #{$testdata}Random.fna RandomN.fna"
runsfxfail "-dna -suf -pl 10 -db #{$testdata}Random.fna"
runsfxfail "-dna -tis -sat plain -db #{$testdata}TTT-small.fna"

Name "gt suffixerator failure"
Keywords "gt_suffixerator"
Test do
  run_test "#{$bin}gt suffixerator -tis -dna -indexname localidx " +
           "-db #{$testdata}Random.fna"
  run_test "#{$bin}gt dev sfxmap -tis -suf -des #{trials()} localidx",
           :retval => 1
end

Name "gt suffixerator bwt"
Keywords "gt_suffixerator"
Test do
  checkbwt(allfiles)
end

allfiles.each do |filename|
  Name "gt suffixerator uint32"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis -indexname sfx -sat uint32 " +
             "-pl -db #{$testdata}#{filename}"
  end
end

1.upto(3) do |parts|
  [0,2].each do |withsmap|
    extra=""
    if withsmap == 1
      extra="-protein"
      extraname="protein"
    elsif withsmap == 2
      extra="-smap TransProt11"
      extraname="TransProt11"
    end
    Name "gt suffixerator+sfxmap protein #{extraname} #{parts} parts"
    Keywords "gt_suffixerator"
    Test do
      checksfx(parts,3,extra,"direct",true,["sw100K1.fsa","sw100K2.fsa"])
    end
  end
end

[true,false].each do |cmp|
  1.upto(2) do |parts|
    ["direct", "bit", "uchar", "ushort", "uint"].each do |sat|
      [0,2].each do |withsmap|
        extra=""
        if withsmap == 1
          extra="-dna"
          extraname="dna"
        elsif withsmap == 2
          extra="-smap TransDNA"
          extraname="TransDNA"
        end
        Name "gt suffixerator+sfxmap dna #{extraname} #{sat} #{parts} parts"
        Keywords "gt_suffixerator"
        Test do
          checksfx(parts,1,extra,sat,cmp,["Random-Small.fna"])
          checksfx(parts,3,extra,sat,cmp,["Random.fna"])
          checksfx(parts,3,extra,sat,cmp,["RandomN.fna"])
          checksfx(parts,2,extra,sat,cmp,["trna_glutamine.fna"])
          checksfx(parts,1,extra,sat,cmp,["TTT-small.fna"])
          checksfx(parts,3,extra,sat,cmp,["RandomN.fna","Random.fna",
                                          "Atinsert.fna"])
        end
      end
    end
  end
end

def checkmapped(args)
  Name "gt suffixerator checkmapped"
  Keywords "gt_suffixerator gttestdata"
  Test do
    run_test "#{$bin}gt suffixerator #{outoptions()} -indexname sfxidx #{args}",
             :maxtime => 600
    run_test "#{$bin}gt dev sfxmap #{outoptions()} #{trials()} -v sfxidx",
             :maxtime => 600
    run_test "#{$bin}gt dev sfxmap #{outoptionsnobck()} -stream -v sfxidx",
             :maxtime => 600
  end
end

def makegreedyfwdmatcall(queryfile,indexarg,ms)
  prog=""
  if ms
    prog="#{$bin}gt matstat -verify"
  else
    prog="#{$bin}gt uniquesub"
  end
  constantargs="-min 1 -max 20 -query #{queryfile} #{indexarg}"
  return "#{prog} -output querypos #{constantargs}"
end

def checkgreedyfwdmat(queryfile,ms)
  run_test makegreedyfwdmatcall(queryfile,"-fmi fmi",ms), :maxtime => 600
  run "mv #{$last_stdout} tmp.fmi"
  run_test makegreedyfwdmatcall(queryfile,"-esa sfx",ms), :maxtime => 600
  run "mv #{$last_stdout} tmp.esa"
  run "diff tmp.esa tmp.fmi"
  run_test makegreedyfwdmatcall(queryfile,"-pck pck",ms), :maxtime => 600
  run "mv #{$last_stdout} tmp.pck"
  run "diff tmp.pck tmp.fmi"
end

def createandcheckgreedyfwdmat(reffile,queryfile)
  run "#{$scriptsdir}/runmkfm.sh #{$bin}/gt 0 . fmi #{reffile}"
  run "#{$bin}gt suffixerator -indexname sfx -tis -suf -dna -v " +
           "-db #{reffile}"
  run "#{$bin}gt packedindex mkindex -tis -indexname pck -db #{reffile} " +
           "-dna -pl -bsize 10 -locfreq 32 -dir rev"
  checkgreedyfwdmat(queryfile,false)
  checkgreedyfwdmat(queryfile,true)
  run "rm -f sfx.* fmi.* pck.*"
end

def grumbach()
  return "#{$gttestdata}DNA-mix/Grumbach.fna/"
end

allfiles.each do |reffile|
  allfiles.each do |queryfile|
    if queryfile != reffile
      Name "gt greedyfwdmat #{reffile} #{queryfile}"
      Keywords "gt_greedyfwdmat small"
      Test do
        createandcheckgreedyfwdmat("#{$testdata}/#{reffile}",
                                   "#{$testdata}/#{queryfile}")
      end 
    end
  end
end

allfiles.each do |reffile|
  Name "gt packedindex #{reffile}"
  Keywords "gt_packedindex small"
  Test do
    run_test "#{$bin}gt packedindex mkindex -tis -indexname pck " +
             " -db #{$testdata}/#{reffile} -dna -pl -bsize 10 " +
             " -locfreq 32 -dir rev", 
             :maxtime => 600
  end
end

if $gttestdata then
  checkmapped("-db " +
              "#{$gttestdata}Iowa/at100K1 " +
              "#{grumbach()}Wildcards.fna " +
              "#{grumbach()}chntxx.fna " +
              "#{grumbach()}hs5hcmvcg.fna " +
              "#{grumbach()}humdystrop.fna " +
              "#{grumbach()}humghcsa.fna " +
              "#{grumbach()}humhdabcd.fna " +
              "#{grumbach()}humhprtb.fna " +
              "#{grumbach()}mipacga.fna " +
              "#{grumbach()}mpocpcg.fna " +
              "#{grumbach()}ychrIII.fna " +
              "-parts 3 -pl")

  checkmapped("-parts 1 -pl -db #{$gttestdata}swissprot/swiss10K " +
              "#{$gttestdata}swissprot/swiss1MB")

  checkmapped("-db #{$gttestdata}swissprot/swiss10K " +
              "#{$gttestdata}swissprot/swiss1MB -parts 3 -pl")

  checkmapped("-parts 2 -pl -smap TransDNA -db  #{$gttestdata}Iowa/at100K1")

  checkmapped("-db #{$gttestdata}swissprot/swiss10K -parts 1 -pl -smap " +
              "TransProt11")
  Name "gt greedyfwdmat at1MB U8"
  Keywords "gt_greedyfwdmat gttestdata"
  Test do
    createandcheckgreedyfwdmat("#{$gttestdata}Iowa/at1MB",
                               "#{$gttestdata}Iowa/U89959.fna")
  end
end
