##' Get protein docking pocket parameters. 
##' Note: Fpocket and MGLTools must be installed and configured correctly. 
##' @title getBoxFromFpocket 
##' @param PDBid PDB ID 
##' @param type FromLigand;FromFpocket;Input 
##' @param PDBboxDb prepared pocket parameter database 
##' @return a dataframe Object 
##' @export 
##' @author Yang Ming 
getPDBbox = function(PDBid,type='FromLigand',center_x=NULL,center_y=NULL,center_z=NULL,size_x=NULL,size_y=NULL,size_z=NULL,python2PTAH='~/mgltools_x86_64Linux2_1.5.6/bin/pythonsh',prepare_ligand4PATH='~/ADtools',PDBboxDb='PDBboxDB.db'){ 
  ##ubuntu 
  ##type:FromLigand;FromFpocket;Input 
  #PDBboxDb first! 
  require(bio3d) 
  require(stringr) 
  tempa=PDBid 
  if (str_detect(tempa,'\\.')){ 
    #tempa=unlist(str_split(tempa,'\\.')) 
    tempa2=basename(normalizePath(tempa)) 
    tempa2=str_replace_all(tempa2,'.pdb','') 
    #tempa2=tempa2[length(tempa2)-1] 
    #tempb=str_sub(tempa,start = nchar(tempa)-3,end = nchar(tempa)) 
    tempPDBid=tolower(tempa2) 
  }else{ 
    tempPDBid=tolower(tempa) 
  } 
  if (type=='FromLigand'){ 
    if (file.exists(PDBboxDb)){ 
      os=paste0('SELECT * FROM PDBboxDB where PDB=','"',tempPDBid,'"') 
      con <- dbConnect(SQLite(),dbname=PDBboxDb) 
      rs <- dbSendQuery(con, os) 
      d1 <- dbFetch(rs) 
      if (nrow(d1)>0){ 
        PDBbox=d1 
      }else{ 
        if (class(try({d2=getBoxFromPdb(pdbFile=tempa,IFprepare=T,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH)},silent =T))!='try-error'){#c(xmean=xmean,ymean=ymean,zmean=zmean,xmax=xmax,ymax=ymax,zmax=zmax) 
          if (sum(!is.na(d2))>0){ 
            d2=as.data.frame(t(d2)) 
            colnames(d2)=c('x','y','z','sizex','sizey','sizez') 
            PDBbox=cbind(PDB=tempPDBid,d2,method='FromLigand') 
          }else{ 
            PDBbox=data.frame() 
          } 
        }else{ 
          PDBbox=data.frame() 
        } 
      } 
    }else{ 
      if (class(try({d2=getBoxFromPdb(pdbFile=tempa,IFprepare=T,python2PTAH=python2PTAH,prepare_ligand4PATH=prepare_ligand4PATH)},silent =T))!='try-error'){#c(xmean=xmean,ymean=ymean,zmean=zmean,xmax=xmax,ymax=ymax,zmax=zmax) 
        if (sum(!is.na(d2))>0){ 
          d2=as.data.frame(t(d2)) 
          colnames(d2)=c('x','y','z','sizex','sizey','sizez') 
          PDBbox=cbind(PDB=tempPDBid,d2,method='FromLigand') 
        }else{ 
          PDBbox=data.frame() 
        } 
      }else{ 
        PDBbox=data.frame() 
      } 
    } 
  }else if (type=='FromFpocket'){ 
    if (file.exists(PDBboxDb)){ 
      os=paste0('SELECT * FROM PDBboxDB where PDB=','"',tempPDBid,'"') 
      con <- dbConnect(SQLite(),dbname=PDBboxDb) 
      rs <- dbSendQuery(con, os) 
      d1 <- dbFetch(rs) 
      if (nrow(d1)>0){ 
        PDBbox=d1 
      }else{ 
        if (class(try({d3=getBoxFromFpocket(pdbID=tempa)},silent = T))!='try-error'){ 
          PDBbox=cbind(PDB=tempPDBid,d3$firstPid[,c('x','y','z','sizex','sizey','sizez')],method='FromFpocket') 
        }else{ 
          PDBbox=data.frame() 
        } 
      } 
    }else{ 
      if (class(try({d3=getBoxFromFpocket(pdbID=tempa)},silent = T))!='try-error'){ 
        PDBbox=cbind(PDB=tempPDBid,d3$firstPid[,c('x','y','z','sizex','sizey','sizez')],method='FromFpocket') 
      }else{ 
        PDBbox=data.frame() 
      } 
    } 
  }else if (type=='Input'){ 
    PDBbox=data.frame(PDB=tempPDBid,x=center_x,y=center_y,z=center_z,sizex=size_x,sizey=size_y,sizez=size_z,method='input') 
  } 
  return(PDBbox) 
} 
 
