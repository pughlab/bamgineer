#from helpers import parameters as params
#
#
#
#def initialize(gaincnv, losscnv):
#    try:
#        
#        if(args.phase):
#            logger.debug(' --- Initialization called --- ')
#            vcf_path = bamgineerHelpers.GetVCF()
#            exons_path = bamgineerHelpers.GetExons()
#            reference_path = bamgineerHelpers.GetRef()
#            
#            
#            vpath, vcf = os.path.split(vcf_path)
#            phasedvcf = "/".join([results_path, sub('.vcf$', '_phased.vcf.gz', vcf)])
#            vcftobed =  "/".join([results_path, sub('.vcf$', '.bed', vcf)])
#            hap1vcf = "/".join([results_path,"hap1_het.vcf"])
#            hap2vcf = "/".join([results_path, "hap2_het.vcf"])
#            hap1vcffiltered = "/".join([results_path, "hap1_het_filtered"])
#            hap2vcffiltered = "/".join([results_path, "hap2_het_filtered"])
#            hap1vcffilteredtobed = "/".join([results_path, "hap1_het_filtered.bed"])
#            hap2vcffilteredtobed = "/".join([results_path, "hap2_het_filtered.bed"])
#            phased_bed =  "/".join([results_path, "PHASED.BED"])
#              
#            phaseVCF(vcf_path, phasedvcf)
#            getVCFHaplotypes(phasedvcf, hap1vcf, hap2vcf)
#            thinVCF(hap1vcf, hap1vcffiltered)
#            thinVCF(hap2vcf, hap2vcffiltered)
#            convertvcftobed(hap1vcffiltered+".recode.vcf", hap1vcffilteredtobed)
#            convertvcftobed(hap2vcffiltered+".recode.vcf", hap2vcffilteredtobed)
#           
#            cmd1 = """sed -i 's/$/\thap1/' """+ hap1vcffilteredtobed
#            cmd2 = """sed -i 's/$/\thap2/' """+ hap2vcffilteredtobed
#            cmd3 = "cat " + hap1vcffilteredtobed + " " + hap2vcffilteredtobed + " > " + 'tmp.bed'
#            cmd4 = "sort -V -k1,1 -k2,2 tmp.bed > " + phased_bed  
#                
#            runCommand(cmd1)
#            runCommand(cmd2)
#            runCommand(cmd3)
#            runCommand(cmd4)
#            os.remove('tmp.bed')  
#            
#            for  event in event_list: 
#                roibed = "/".join([haplotype_path,  event + "_roi.bed"])
#                exonsinroibed = "/".join([haplotype_path,   event + "_exons_in_roi.bed"])
#                exonsinhetbed = "/".join([haplotype_path,  event + "_exons_in_het.bed"])
#                nonhetbed = "/".join([haplotype_path, event + "_non_het.bed"])
#                hetbed = "/".join([haplotype_path, event + "_het.bed"])
#                hetsnpbed = "/".join([haplotype_path,  event + "_het_snp.bed"])
#                
#                logger.debug("Calling initialization for " + str(event))
#        
#                #intersectBed( exons_path, locals()[event + 'cnv'], exonsinroibed, wa=True)
#                #intersectBed( exonsinroibed, phased_bed, exonsinhetbed, wa=True)  
#                #intersectBed(phased_bed, exonsinroibed, hetsnpbed, wa=True)
#                # 
#                #subtractBeds(exonsinroibed, phased_bed , nonhetbed) # so that 
#                #intersectBed( phased_bed, exonsinroibed, hetbed, wa=True)
#                #
#                #splitBed(hetsnpbed, event+'_het_snp_')
#                #splitBed(hetbed, event+'_het_')
#                #splitBed(nonhetbed, event+'_non_het_')             
#            
#    except:
#        logger.error('Initialization error ', sys.exc_info()[0])
#        raise
#    return
