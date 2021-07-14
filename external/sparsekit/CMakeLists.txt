##########################################################################
# sparsekit directory
##########################################################################
set(files
        BLASSM/blassm.f 
        FORMATS/chkfmt1.f       
        FORMATS/chkun.f         
        FORMATS/formats.f       
        FORMATS/rvbr.f          
        FORMATS/unary.f         
        INOUT/chkio.f      
        INOUT/hb2pic.f     
        INOUT/hb2ps.f      
        INOUT/inout.f      
        ORDERINGS/ccn.f    
        ORDERINGS/color.f  
        ORDERINGS/dsepart.f
        ITSOL/ilut.f        
        ITSOL/itaux.f       
        ITSOL/iters.f       
        ITSOL/rilut.f       
        ITSOL/riter2.f      
        ITSOL/riters.f)     

add_library(sll_sparsekit ${files})
##########################################################################
