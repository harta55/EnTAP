/*
 *
 * Developed by Alexander Hart
 * Plant Computational Genomics Lab
 * University of Connecticut
 *
 * For information, contact Alexander Hart at:
 *     entap.dev@gmail.com
 *
 * Copyright 2017-2018, Alexander Hart, Dr. Jill Wegrzyn
 *
 * This file is part of EnTAP.
 *
 * EnTAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * EnTAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with EnTAP.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ENTAP_EGGNOGLEVELS_H
#define ENTAP_EGGNOGLEVELS_H

#include <map>

std::map<std::string,std::string> EGGNOG_LEVELS = {
        {"acoNOG", "Aconoidasida"},
	    {"agaNOG", "Agaricales"},
	    {"agarNOG", "Agaricomycetes"},
	    {"meNOG", "Animals"},
	    {"apiNOG", "Apicomplexa"},
	    {"arthNOG", "Arthrodermataceae"},
	    {"artNOG", "Arthropoda"},
        {"ascNOG", "Ascomycota"},
        {"aveNOG", "Aves"},
        {"basNOG", "Basidiomycota"},
        {"biNOG", "Bilateria"},
        {"braNOG", "Brassicales"},
        {"carNOG", "Carnivora"},
        {"chaNOG", "Chaetomiaceae"},
        {"chloroNOG", "Chlorophyta"},
        {"chorNOG", "Chordata"},
        {"chrNOG", "Chromadorea"},
        {"cocNOG", "Coccidia"},
        {"cryNOG", "Cryptosporidiidae"},
        {"debNOG", "Debaryomycetaceae"},
        {"dipNOG", "Diptera"},
        {"dotNOG", "Dothideomycetes"},
        {"droNOG", "Drosophilidae"},
        {"euNOG", "Eukaryotes"},
        {"eurotNOG", "Eurotiales"},
        {"euroNOG", "Eurotiomycetes"},
        {"fiNOG", "Fishes"},
        {"fuNOG", "Fungi"},
        {"haeNOG", "Haemosporida"},
        {"homNOG", "Hominidae"},
        {"hymNOG", "Hymenoptera"},
        {"hypNOG", "Hypocreales"},
        {"inNOG", "Insects"},
        {"kinNOG", "Kinetoplastida"},
        {"lepNOG", "Lepidoptera"},
        {"lilNOG", "Liliopsida"},
        {"magNOG", "Magnaporthales"},
        {"maNOG", "Mammals"},
        {"necNOG", "Nectriaceae"},
        {"nemNOG", "Nematodes"},
        {"onyNOG", "Onygenales"},
        {"opiNOG", "Opisthokonts"},
        {"perNOG", "Peronosporales"},
        {"pleNOG", "Pleosporales"},
        {"poaNOG", "Poales"},
        {"prNOG", "Primates"},
        {"rhaNOG", "Rhabditida"},
        {"roNOG", "Rodents"},
        {"sacNOG", "Saccharomycetaceae"},
        {"saccNOG", "Saccharomycetes"},
        {"sorNOG", "Sordariales"},
        {"sordNOG", "Sordariomycetes"},
        {"strNOG", "Streptophyta"},
        {"spriNOG", "Supraprimates"},
        {"treNOG", "Tremellales"},
        {"veNOG", "Vertebrates"},
        {"virNOG", "Viridiplantae"},
        {"bactNOG", "Bacteria"},
        {"bctoNOG", "Bacteroidetes"},
        {"chlNOG", "Chlorobi"},
        {"cyaNOG", "Cyanobacteria"},
        {"proNOG", "Proteobacteria"},
        {"gproNOG", "Proteobacteria_gamma"},
        {"firmNOG", "Firmicutes"},
        {"deiNOG", "Deinococcusthermus"},
        {"aproNOG", "Proteobacteria_alpha"},
        {"bproNOG", "Proteobacteria_beta"},
        {"dproNOG", "Proteobacteria_delta"},
        {"eproNOG", "Proteobacteria_epsilon"},
        {"chlorNOG", "Chloroflexi"},
        {"fusoNOG", "Fusobacteria"},
        {"aciNOG", "Acidobacteria"},
        {"delNOG", "delta/epsilon"},
        {"verNOG", "Verrucomicrobia"},
        {"bacNOG", "Bacilli"},
        {"flaNOG", "Flavobacteriia"},
        {"sphNOG", "Sphingobacteriia"},
        {"cloNOG", "Clostridia"},
        {"bacteNOG", "Bacteroidia"},
        {"aquNOG", "Aquificae"},
        {"chloNOG", "Chloroflexi"},
        {"therNOG", "Thermotogae"},
        {"defNOG", "Deferribacteres"},
        {"actNOG", "Actinobacteria"},
        {"verrNOG", "Verrucomicrobiae"},
        {"plaNOG", "Planctomycetes"},
        {"spiNOG", "Spirochaetes"},
        {"chlaNOG", "Chlamydiae"},
        {"acidNOG", "Acidobacteriia"},
        {"dehNOG", "Dehalococcoidetes"},
        {"synNOG", "Synergistetes"},
        {"eryNOG", "Erysipelotrichi"},
        {"tenNOG", "Tenericutes"},
        {"cytNOG", "Cytophagia"},
        {"negNOG", "Negativicutes"},
        {"arNOG", "Archaea"},
        {"creNOG", "Crenarchaeota"},
        {"eurNOG", "Euryarchaeota"},
        {"metNOG", "Methanobacteria"},
        {"methNOG", "Methanococci"},
        {"halNOG", "Halobacteria"},
        {"theNOG", "Thermoplasmata"},
        {"thermNOG", "Thermococci"},
        {"arcNOG", "Archaeoglobi"},
        {"methaNOG", "Methanomicrobia"},
        {"thaNOG", "Thaumarchaeota"},
        {"NOG", "Ancestor"}
};


#endif //ENTAP_EGGNOGLEVELS_H
