setwd('C:/Users/jdx66/Box Sync/job/incubator/county_level_estimates')
library(sp);
library(rgeos)
library(maptools);
library(RColorBrewer)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(ggmap)
library(stringr)
library(sas7bdat)
library(Rcpp)
library(MCMCpack)
library(HI)


map_mo = map_data('county','missouri')
getLabelPoint <- # Returns a county-named list of label points
  function(county) {Polygon(county[c('long', 'lat')])@labpt}
centroids <- by(map_mo, map_mo$subregion, getLabelPoint)
centroids <- do.call("rbind.data.frame", centroids)

names(centroids) <- c('long', 'lat')

countyName = read.table('data/countyName.txt',sep=';')
countyName[,1] = as.character(countyName[,1])
get_county_fips = function(x){
  temp = unlist(strsplit(x,'='))
  temp = substr(temp[2], 3, 5)
  return(as.numeric(temp))
}
get_county_fips(countyName[2,1])
fips = as.numeric(sapply(countyName[,1],get_county_fips))
names = tolower(countyName[,2])

get_fips_data = function(map_mo,countyName,fips){
  n = length(countyName)
  df = data.frame(orders=1:n,
                  fips = fips,
                  subregion=countyName,
                  label = paste(1:n,countyName,sep='\n '))
  fips_data = merge(map_mo, df, sort = FALSE, by = "subregion")
  fips_data = fips_data[order(fips_data$order),]
  return(fips_data)
}

get_centroids_fips_df = function(centroids_df,countyName, fips){
  temp = centroids_df[countyName,]
  n = length(countyName)
  temp$label = paste(1:n,countyName,fips,sep='\n ')
  return(temp)
}
df_fips = get_fips_data(map_mo,names, fips)
df_fips_text = get_centroids_fips_df(centroids,names, fips)

# BRFSS regions and counties
region2county = list()
#1. Kansas city Metro
region2county[[1]] = c(165,49,25,47,177,95,107,37,13)
#2. St.Louis Metro
region2county[[2]] = c(113,219,183,71,189,510,99,221)
#3.Central
region2county[[3]] = c(159,89,53,141,29,105,19,135,131,169,51,7,27,151,125,161,139,73,55,65)
#4.Southwest
region2county[[4]] = c(217,11,97,145,119,83,185,39,57,109,9,15,85,167,77,43,213,59,225,209)
#5.Southest
region2county[[5]] = c(229,67,153,215,91,203,149,93,179,35,181,186,187,123,223,23,157,17,31,201,207,133,143,69,155)
#6. Northwest
region2county[[6]] = c(5,87,147,3,21,227,75,63,81,61,33,101)
#7. Northeast
region2county[[7]] = c(129,79,117,171,211,115,41,195,197,1,121,175,199,103,205,137,45,111,127,173,163)

#note, since fips is not countiuous, the following list will have NULLs
county2region = list()
county2region_df = data.frame(fips=fips, brfss_region=0)
for(i in 1:length(region2county)){
  region = i;
  counties = region2county[[region]]
  n_counties = length(counties)
  for(j in 1:n_counties){
    county = counties[j]
    county2region[county] = region
    
    county2region_df[county2region_df$fips==county,2] = region
  }
}

n <- 7
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(21)
col_vector_sample = sample(col_vector, n)

df_fips_regions = merge(x=df_fips, y=county2region_df, by=c('fips'))
brfss_fips_plot = ggplot(df_fips_regions,aes(long,lat))+
  geom_polygon(aes(group=group,fill=as.factor(brfss_region)),colour='black')+
  scale_fill_manual(values=col_vector_sample,
                    guide = guide_legend(title = 'BRFSS\nregions'))+
  geom_text(data=df_fips_text, aes(x=long,y=lat,label=label),size=1.5)+
  ggtitle('MO FIPS MAP')+
  theme_bw()+
  theme(plot.title = element_text(size=12),
        strip.text.x = element_text(size = 8, colour = "red", angle = 0))+
  labs(x="",y="")+
  coord_map()

brfss2012 = read.sas7bdat(file='data/saebrfss2012.sas7bdat')

sub_cols = c('CTYCODE1', 'X_IMPAGE', 'SEX', 'X_IMPRACE', 'HADSIGM3', 'X_RFSIGM2')
sub_brfss = brfss2012[,sub_cols]

colnames(sub_brfss) = c('county', 'age', 'sex', 'race', 'hs', 'hs50')

#sex:1=male, 2=female
#generate age group
n_obs = dim(sub_brfss)[1]
age_group = c()
race_group = c()
region = c()
class_vec = c()
for(i in 1:n_obs){
  cty = sub_brfss[i,'county']
  age = sub_brfss[i, 'age']
  race = sub_brfss[i, 'race']
  sex_i = sub_brfss[i, 'sex']
  if(cty%in%fips){
    region_i = county2region[[cty]]
  }else{
    region_i=NA
  }
  
  if(age<50 & age>=18){
    age_group_i = 0
  }else if(age>=50 & age<65){
    age_group_i = 1
  }else if (age>=65 & age<75){
    age_group_i = 2
  }else if (age>=75 & age<200){
    age_group_i = 3
  }else{
    age_group_i = 3
  }
  
  if(race==1){
    race_group_i = 1
  }else{
    race_group_i = 2
  }
  
  class_i = paste(sex_i, age_group_i, race_group_i, sep='_')
  region = c(region, region_i)
  race_group = c(race_group, race_group_i)
  age_group = c(age_group, age_group_i)
  class_vec = c(class_vec, class_i)
}

sub_brfss$'region'=region
sub_brfss$'race_group'=race_group
sub_brfss$'age_group'=age_group
sub_brfss$'class_vec'=class_vec

brfss_old = sub_brfss[sub_brfss$age_group!=0,]

class_info = table(brfss_old$class_vec)

#"1_1_1" "1_1_2" "1_2_1" "1_2_2" "1_3_1" "1_3_2" "2_1_1" "2_1_2" "2_2_1" "2_2_2" "2_3_1" "2_3_2"
class_names = names(class_info)


#remove unknown counties


unknown_county = brfss_old$county>550 | is.na(brfss_old$county)

brfss_old = brfss_old[!unknown_county,]

#check hs, hs50 condition

#remove unknown hs50 since it helped nothing with our model
brfss_old = brfss_old[!is.na(brfss_old$hs50),]
#recode hs50 to be 0(no screen), 1(has screened)
brfss_old$'hs50ind' = as.numeric(brfss_old$hs50==1)
#[1] 3807   10
#read population data
population = read.csv(file='data/population_by_class.txt',header = FALSE)
colnames(population) = c('county', 'sex', 'age_group', 'race_group', 'rate', 'count', 'pop')
#sex: 0 is male, 1 is female
#age: 0=50-64, 1=65-74, 2=75+
#race:0=white, 1=non-white
#county:0-114
#recode the variables to match prior definations
population$county_re = population$county+1 #1-115
population$sex_re = population$sex + 1
population$age_re = population$age_group + 1
population$race_re = population$race_group + 1
county_fips = numeric(dim(population)[1])
for(i in 1:nrow(population)){
  county_fips[i] = fips[population$county_re[i]]
}
population$'county_fips' = county_fips
class_vec_2 = numeric(dim(population)[1])
for(i in 1:nrow(population)){
  class_vec_2[i] = paste(population$sex_re[i], population$age_re[i], population$race_re[i], sep='_')
}
population$'class_vec' = class_vec_2
#generate yij, nij, Nij mat, and R prior
y_mat = n_mat = N_mat = matrix(NA, length(fips), length(class_names))
head(brfss_old)
table(brfss_old$county)
for(i in 1:length(fips)){
  cty = fips[i]
  sub_data = brfss_old[brfss_old$county==cty,]
  sub_pop = population[population$county_fips==cty, ]
  for(j in 1:length(class_names)){
    class_id = class_names[j]
    subsub_data = sub_data[sub_data$class_vec==class_id,]
    subsub_pop = sub_pop[sub_pop$class_vec==class_id,]
    
    if(nrow(subsub_pop)==1){
      N_ij = as.numeric(subsub_pop$pop)
    }else{
      print('Invalid N_ij')
    }
    
    if(dim(subsub_data)[1]>=1){
      n_ij = nrow(subsub_data)
      y_ij = sum(subsub_data$hs50ind)
    }else{
      #cat(i,j,'\n',sep='-')
      n_ij = 0
      y_ij = 0
    }
    
    if(N_ij==0){
      y_ij = 0
    }
    if(N_ij<n_ij){
      y_ij = 0
    }
    
    #added, since binomial regression can't run with n=0, make it n=1
    N_ij_modifiy = N_ij
    n_ij_modifiy = n_ij
    #     if(N_ij==0){
    #       N_ij_modifiy=1
    #     }
    #     if(n_ij==0){
    #       #cat(i,j,'\n',sep='-')
    #       n_ij_modifiy=1
    #     }
    
    
    y_mat[i, j] = y_ij
    n_mat[i, j] = n_ij_modifiy
    N_mat[i, j] = N_ij_modifiy
  }
}
#construct ragged data
head(brfss_old)
dim(brfss_old)
class_names
n_class = length(class_names)
n_county = length(fips)
n_data = nrow(brfss_old)
class_to_label = list()
for(i in 1:n_class){
  class_to_label[[class_names[i]]] = i
}
#class_i = paste(sex_i, age_group_i, race_group_i, sep='_')
sex_vec = race_vec = age_vec = numeric(n_class)
for(i in 1:n_class){
  the_name = class_names[i]
  temp = as.numeric(strsplit(the_name, '_')[[1]])
  sex_vec[i] = temp[1]
  age_vec[i] = temp[2]
  race_vec[i] = temp[3]
}#for winbugs use
#generate a dic
class_dic = data.frame(label = 1:length(class_names),class_names=class_names, sex_vec, age_vec, race_vec)


county_order = numeric(n_data)
class_label = numeric(n_data)
for(i in 1:n_data){
  cty = brfss_old$county[i]
  county_order[i] = which(fips==cty)
  class_label[i] = class_to_label[[brfss_old$class_vec[i]]]
}
brfss_old$'county_order' = county_order
brfss_old$'class_label' = class_label
#construc ragged data
data_list = list()
counties_no_info = c()
for(i in 1:n_county){
  data_i = brfss_old[brfss_old$county_order==i,]
  n_i = nrow(data_i)
  if(n_i>0){
    temp1 = aggregate(data_i$hs50ind, by=list(class_label=data_i$class_label), FUN=sum)
    temp2 = aggregate(rep(1,n_i), by=list(class_label=data_i$class_label), FUN=sum)
    temp = merge(temp1,temp2,by='class_label')
    colnames(temp) = c('class_label','y','n')
    temp$county = i
    region_i = data_i$region
    if(sum(region_i-min(region_i))!=0){
      print('error region!')
    }else{
      temp$region = min(region_i)
    }
    data_list[[i]] = temp
  }else{
    counties_no_info = c(counties_no_info, i)
  }
}

data_ragged = do.call('rbind', data_list)

Nrest_mat = N_mat - n_mat
N_county = apply(N_mat, 1, sum)
y_insample = apply(y_mat, 1, sum)
#county attributes
county_attributes = read.table(file='data/county_attributes.txt',header = FALSE, sep='\t')
county_attributes = county_attributes[1:115,]
colnames(county_attributes) = c('county_raw', 'below_poperty', 'below_high_school',
                                'median_income','crc_sae','above_bach', 'below_9th')

#Reformat
county_attributes$below_poperty = as.numeric(as.character(county_attributes$below_poperty))/100
county_attributes$below_high_school = as.numeric(as.character(county_attributes$below_high_school))/100
county_attributes$median_income = as.numeric(as.character(county_attributes$median_income))/100#to thousands dollors
county_attributes$crc_sae = as.numeric(as.character(county_attributes$crc_sae))/100
county_attributes$above_bach = as.numeric(as.character(county_attributes$above_bach))/100
county_attributes$below_9th = as.numeric(as.character(county_attributes$below_9th))/100
county_attributes$county = 1:115
county_attributes$county_raw = as.character(county_attributes$county_raw)
head(county_attributes)
#make sure county is in the right order
get_county_fips_2 = function(x){
  temp = unlist(strsplit(x,';'))
  return(as.numeric(temp[1]))
}
#cut continous to category for possible future use
cut_cont = function(x){
  cut_points = quantile(x, probs=c(0, 0.2, 0.4, 0.6, 0.8, 1))
  cut_points[1] = 0.95*cut_points[1]
  bins = cut(x, cut_points, labels=c(1:5))
  result = as.numeric(as.character(bins))
  return(result)
}
county_attributes$below_poperty_cat = cut_cont(county_attributes$below_poperty)
county_attributes$below_high_school_cat = cut_cont(county_attributes$below_high_school)
county_attributes$median_income_cat = cut_cont(county_attributes$median_income)
county_attributes$above_bach_cat = cut_cont(county_attributes$above_bach)
county_attributes$below_9th_cat = cut_cont(county_attributes$below_9th)
head(county_attributes)
#one is ragged, one is not;
data_ragged_contAttr = merge(x=data_ragged, y=county_attributes, by='county', all.x=TRUE)
brfss_old_contATTR = merge(x=brfss_old, y=county_attributes, by.x='county_order', by.y='county', all.x=TRUE)

#adjacency matrix
#try agg
num <- c(7, 5, 2, 7, 4, 4, 4, 6, 5, 7, 
         4, 5, 6, 5, 7, 4, 6, 6, 4, 5, 
         7, 6, 3, 4, 6, 5, 6, 6, 6, 6, 
         6, 6, 6, 7, 4, 7, 7, 6, 6, 6, 
         5, 6, 5, 3, 5, 5, 7, 5, 4, 5, 
         5, 6, 6, 6, 7, 4, 4, 6, 6, 2, 
         7, 5, 5, 4, 4, 7, 2, 5, 6, 7, 
         6, 5, 4, 5, 4, 6, 3, 2, 5, 7, 
         6, 4, 3, 6, 6, 4, 4, 6, 6, 5, 
         3, 4, 7, 3, 6, 4, 6, 3, 4, 4, 
         6, 5, 7, 4, 5, 4, 8, 4, 5, 5, 
         7, 6, 3, 4, 1)
adj <- c(
  105, 99, 98, 86, 61, 58, 52, 
  74, 44, 38, 32, 11, 
  74, 44, 
  88, 87, 82, 70, 69, 14, 10, 
  104, 73, 60, 55, 
  108, 49, 29, 20, 
  108, 93, 42, 19, 
  93, 80, 71, 43, 42, 15, 
  111, 103, 79, 62, 16, 
  88, 68, 45, 27, 26, 14, 4, 
  83, 32, 25, 2, 
  111, 103, 91, 35, 18, 
  89, 59, 32, 31, 25, 17, 
  76, 70, 26, 10, 4, 
  85, 71, 66, 53, 43, 30, 8, 
  103, 100, 79, 9, 
  97, 89, 59, 54, 21, 13, 
  111, 101, 91, 90, 75, 12, 
  51, 48, 42, 7, 
  108, 93, 84, 29, 6, 
  97, 88, 61, 59, 58, 45, 17, 
  112, 106, 104, 55, 39, 34, 
  99, 56, 52, 
  89, 83, 48, 25, 
  89, 83, 32, 24, 13, 11, 
  76, 68, 66, 14, 10, 
  97, 80, 71, 68, 45, 10, 
  110, 81, 47, 37, 36, 33, 
  84, 55, 49, 39, 20, 6, 
  112, 84, 53, 43, 39, 15, 
  59, 41, 40, 38, 32, 13, 
  38, 31, 25, 13, 11, 2, 
  107, 101, 90, 81, 47, 28, 
  114, 112, 107, 106, 77, 46, 22, 
  103, 78, 72, 12, 
  110, 109, 96, 92, 50, 37, 28, 
  109, 81, 76, 70, 63, 36, 28, 
  113, 74, 41, 32, 31, 2, 
  112, 84, 55, 30, 29, 22, 
  105, 65, 59, 58, 41, 31, 
  113, 65, 40, 38, 31, 
  93, 80, 51, 19, 8, 7, 
  93, 84, 30, 15, 8, 
  74, 3, 2, 
  97, 88, 27, 21, 10, 
  107, 101, 77, 75, 34, 
  111, 110, 95, 90, 62, 33, 28, 
  89, 54, 51, 24, 19, 
  73, 55, 29, 6, 
  110, 96, 95, 94, 36, 
  80, 54, 48, 42, 19, 
  102, 99, 61, 56, 23, 1, 
  114, 112, 107, 85, 30, 15, 
  97, 89, 80, 51, 48, 17, 
  104, 73, 49, 39, 29, 22, 5, 
  102, 64, 52, 23, 
  109, 92, 82, 70, 
  105, 61, 59, 40, 21, 1, 
  58, 40, 31, 21, 17, 13, 
  73, 5, 
  102, 88, 69, 58, 52, 21, 1, 
  111, 95, 79, 47, 9, 
  85, 81, 76, 66, 37, 
  102, 87, 69, 56, 
  105, 86, 41, 40, 
  85, 76, 71, 68, 63, 26, 15, 
  100, 72, 
  71, 66, 27, 26, 10, 
  102, 88, 87, 64, 61, 4, 
  109, 82, 76, 57, 37, 14, 4, 
  80, 68, 66, 27, 15, 8, 
  103, 100, 78, 67, 35, 
  60, 55, 49, 5, 
  113, 44, 38, 3, 2, 
  101, 91, 46, 18, 
  70, 66, 63, 37, 26, 14, 
  106, 46, 34, 
  72, 35, 
  95, 94, 62, 16, 9, 
  97, 71, 54, 51, 42, 27, 8, 
  107, 85, 63, 37, 33, 28, 
  87, 70, 57, 4, 
  25, 24, 11, 
  93, 43, 39, 30, 29, 20, 
  107, 81, 66, 63, 53, 15, 
  105, 98, 65, 1, 
  82, 69, 64, 4, 
  69, 61, 45, 21, 10, 4, 
  54, 48, 25, 24, 17, 13, 
  111, 101, 47, 33, 18, 
  75, 18, 12, 
  109, 96, 57, 36, 
  108, 84, 43, 42, 20, 8, 7, 
  95, 79, 50, 
  110, 94, 79, 62, 50, 47, 
  115, 92, 50, 36, 
  80, 54, 45, 27, 21, 17, 
  99, 86, 1, 
  98, 52, 23, 1, 
  103, 72, 67, 16, 
  107, 90, 75, 46, 33, 18, 
  69, 64, 61, 56, 52, 
  111, 100, 72, 35, 16, 12, 9, 
  106, 55, 22, 5, 
  86, 65, 58, 40, 1, 
  104, 77, 34, 22, 
  114, 101, 85, 81, 53, 46, 34, 33, 
  93, 20, 7, 6, 
  92, 70, 57, 37, 36, 
  95, 50, 47, 36, 28, 
  103, 90, 62, 47, 18, 12, 9, 
  114, 53, 39, 34, 30, 22, 
  74, 41, 38, 
  112, 107, 53, 34, 
  96);
sumNumNeigh <- 588;
W       <- matrix(0, 115,115);
idx     <- 1;
neighs  <- 1;
for (i in 1:sumNumNeigh) {
  W[idx, adj[i]] <- 1;
  W[adj[i], idx] <- 1;
  if (neighs < num[idx]) {
    neighs <- neighs + 1;
  } else {
    idx     <- idx + 1;
    neighs  <- 1;
  }
}


















































