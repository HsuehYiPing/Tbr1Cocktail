import csv
import os

dsets = ["NO.1", "NO.2", "NO.3", "NO.8", "NO.13", "NO.15", "NO.16", "NO.20"]
dtreatments = ["Water", "Cocktail"]
dstatus = ["OE", "RSI"]

for dest in dsets:
    # read Multi-session_registration file
    cellOE = {}
    cellRSI = {}
    cellOE2 = {}
    cellRSI2 = {}
    cells = []
    with open('input_files/Multi-session_registration_' + dest + '.csv', newline='') as regfile:
        rows = csv.DictReader(regfile)
        for row in rows:
            cells.append(row['label'])
            if row['Water_OE'] is not None: cellOE[row['Water_OE']] = row['label']
            if row['Water_RSI'] is not None: cellRSI[row['Water_RSI']] = row['label']
            if row['Cocktail_OE'] is not None: cellOE2[row['Cocktail_OE']] = row['label']
            if row['Cocktail_RSI'] is not None: cellRSI2[row['Cocktail_RSI']] = row['label']
    regfile.close()

    # read each output file
    val = {}
    for dt in dtreatments:
        val[dt]={}
        for st in dstatus:
            val[dt][st]={}
            resfile = "input_files/" + dt + "/" + dest + "/" + st + "/output_files/summary_P.csv"
            with open(resfile, newline='') as infile:
                rows = csv.DictReader(infile)
                for row in rows:
                    if (row['behavior correlation *'] == "inhibition" or row['behavior correlation *'] == "activation"):
                        if dt == "Water":
                            if st == "OE": 
                                val[dt][st][cellOE[row['cell ID']]] = row['behavior correlation *']
                            elif st == "RSI": 
                                val[dt][st][cellRSI[row['cell ID']]] = row['behavior correlation *']
                        if dt == "Cocktail":
                            if st == "OE": 
                                val[dt][st][cellOE2[row['cell ID']]] = row['behavior correlation *']
                            elif st == "RSI": 
                                val[dt][st][cellRSI2[row['cell ID']]] = row['behavior correlation *']
            infile.close()
            print(dt + "_" + st + "_" + dest + " done")

    # output cell table
    outfile1 = "parsing_cell_type_to_links_" + dest + ".txt"
    outfile2 = "parsing_cell_type_to_links_" + dest + "_nodes.txt"

    with open(outfile1, 'w') as fout1, open(outfile2, 'w') as fout2:
        # header
        colname = []
        for dt in dtreatments:
            for st in dstatus:
                sname = dt + "_" + st
                colname.append(sname)
        cname = "cell" + "\t" + '\t'.join(colname) + "\n"
        fout1.write(cname)

        # row
        for c in cells:
            rowvals = []
            for dt in dtreatments:
                for st in dstatus:
                    if c in val[dt][st]:
                        rowvals.append(val[dt][st][c])
                        rowval2 = dt + "\t" + st + "\t" + c + "\t" + val[dt][st][c] + "\n"
                        fout2.write(rowval2)
                    else:
                        rowvals.append("")
            rowval = c + "\t" + '\t'.join(rowvals) + "\n"
            fout1.write(rowval)
            
            
    ######################################################################
    # output similarity table
    # read each output file
    val2 = {}
    for dt in dtreatments:
        val2[dt]={}
        for st in dstatus:
            val2[dt][st]={}
            resfile = "input_files/" + dt + "/" + dest + "/" + st + "/output_files/summary_P.csv"
            with open(resfile, newline='') as infile:
                rows = csv.DictReader(infile)
                for row in rows:
                    if dt == "Water":
                        if st == "OE": 
                            val2[dt][st][cellOE[row['cell ID']]] = row[list(row.keys())[1]]
                        elif st == "RSI": 
                            val2[dt][st][cellRSI[row['cell ID']]] = row[list(row.keys())[1]]
                    if dt == "Cocktail":
                        if st == "OE": 
                            val2[dt][st][cellOE2[row['cell ID']]] = row[list(row.keys())[1]]
                        elif st == "RSI": 
                            val2[dt][st][cellRSI2[row['cell ID']]] = row[list(row.keys())[1]]
            infile.close()
            print(dt + "_" + st + "_" + dest + "_similarity done")
            
            
    # output cell table
    outfile3 = "parsing_cell_type_to_links_" + dest + "_similarity.txt"
    outfile4 = "parsing_cell_type_to_links_" + dest + "_nodes_similarity.txt"

    with open(outfile3, 'w') as fout3, open(outfile4, 'w') as fout4:
        # header
        colname = []
        for dt in dtreatments:
            for st in dstatus:
                sname = dt + "_" + st
                colname.append(sname)
        cname = "cell" + "\t" + '\t'.join(colname) + "\n"
        fout3.write(cname)

        # row
        for c in cells:
            rowvals = []
            for dt in dtreatments:
                for st in dstatus:
                    if c in val2[dt][st]:
                        rowvals.append(val2[dt][st][c])
                        rowval2 = dt + "\t" + st + "\t" + c + "\t" + val2[dt][st][c]+ "\n"
                        fout4.write(rowval2)
                    else:
                        rowvals.append("")
            rowval = c + "\t" + '\t'.join(rowvals) + "\n"
            fout3.write(rowval)
            
            
    ######################################################################
    # output experience dependent cell table
    # read each output file
    val3 = {}
    val4 = {}
    for dt in dtreatments:
        val3[dt]={}
        val4[dt]={}
        for st in dstatus:
            val3[dt][st]={}
            val4[dt][st]={}
            resfile = "input_files/" + dt + "/" + dest + "/" + st + "/output_files/summary_P.csv"
            with open(resfile, newline='') as infile:
                rows = csv.DictReader(infile)
                for row in rows:
                    if dt == "Water":
                        if st == "OE": 
                            val3[dt][st][cellOE[row['cell ID']]] = row["pearson_r"]
                            val4[dt][st][cellOE[row['cell ID']]] = row["Experience_dependent_cell_type"]
                        elif st == "RSI": 
                            val3[dt][st][cellRSI[row['cell ID']]] = row["pearson_r"]
                            val4[dt][st][cellRSI[row['cell ID']]] = row["Experience_dependent_cell_type"]
                    if dt == "Cocktail":
                        if st == "OE": 
                            val3[dt][st][cellOE2[row['cell ID']]] = row["pearson_r"]
                            val4[dt][st][cellOE2[row['cell ID']]] = row["Experience_dependent_cell_type"]
                        elif st == "RSI": 
                            val3[dt][st][cellRSI2[row['cell ID']]] = row["pearson_r"]
                            val4[dt][st][cellRSI2[row['cell ID']]] = row["Experience_dependent_cell_type"]
            infile.close()
            print(dt + "_" + st + "_" + dest + "_Experience_dependent_cell done")
            
            
        # output cell table
        outfile5 ="parsing_cell_type_to_links_" + dest + "_pearson_r.txt"
        outfile6 = "parsing_cell_type_to_links_" + dest + "_nodes_pearson_r.txt"
        outfile7 = "parsing_cell_type_to_links_" + dest + "_experience_dependent_cells.txt"

    with open(outfile5, 'w') as fout5, open(outfile6, 'w') as fout6, open(outfile7, 'w') as fout7:
        # header
        colname = []
        for dt in dtreatments:
            for st in dstatus:
                sname = dt + "_" + st
                colname.append(sname)
        cname = "cell" + "\t" + '\t'.join(colname) + "\n"
        fout5.write(cname)
        fout7.write(cname)

        # row
        for c in cells:
            rowval3s = []
            rowval4s = []
            for dt in dtreatments:
                for st in dstatus:
                    if c in val3[dt][st]:
                        rowval3s.append(val3[dt][st][c])
                        rowval2 = dt + "\t" + st + "\t" + c + "\t" + val3[dt][st][c] + "\n"
                        fout6.write(rowval2)
                        rowval4s.append(val4[dt][st][c])
                    else:
                        rowval3s.append("")
                        rowval4s.append("")
            rowval3 = c + "\t" + '\t'.join(rowval3s) + "\n"
            rowval4 = c + "\t" + '\t'.join(rowval4s) + "\n"
            fout5.write(rowval3)
            fout7.write(rowval4)
            
            
    ######################################################################
    # output mean activity table
    # read each output file
    val5 = {}
    for dt in dtreatments:
        val5[dt]={}
        for st in dstatus:
            val5[dt][st]={}
            resfile = "input_files/" + dt + "/" + dest + "/" + st + "/output_files/summary_P.csv"
            with open(resfile, newline='') as infile:
                rows = csv.DictReader(infile)
                for row in rows:
                    if dt == "Water":
                        if st == "OE": 
                            val5[dt][st][cellOE[row['cell ID']]] = row["behavior_activity_mean"]
                        elif st == "RSI": 
                            val5[dt][st][cellRSI[row['cell ID']]] = row["behavior_activity_mean"]
                    if dt == "Cocktail":
                        if st == "OE": 
                            val5[dt][st][cellOE2[row['cell ID']]] = row["behavior_activity_mean"]
                        elif st == "RSI": 
                            val5[dt][st][cellRSI2[row['cell ID']]] = row["behavior_activity_mean"]
            infile.close()
            print(dt + "_" + st + "_" + dest + "_mean_activity done")
            
            
    # output cell table
    outfile8 = "parsing_cell_type_to_links_" + dest + "_mean.txt"
    outfile9 = "parsing_cell_type_to_links_" + dest + "_nodes_mean.txt"

    with open(outfile8, 'w') as fout8, open(outfile9, 'w') as fout9:
        # header
        colname = []
        for dt in dtreatments:
            for st in dstatus:
                sname = dt + "_" + st
                colname.append(sname)
        cname = "cell" + "\t" + '\t'.join(colname) + "\n"
        fout8.write(cname)

        # row
        for c in cells:
            rowvals = []
            for dt in dtreatments:
                for st in dstatus:
                    if c in val5[dt][st]:
                        rowvals.append(val5[dt][st][c])
                        rowval2 = dt + "\t" + st + "\t" + c + "\t" + val5[dt][st][c]+ "\n"
                        fout9.write(rowval2)
                    else:
                        rowvals.append("")
            rowval = c + "\t" + '\t'.join(rowvals) + "\n"
            fout8.write(rowval)