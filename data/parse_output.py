import sys


def main():
    fn=sys.argv[1]

    file=open(fn)
    num_infer={}
    r_infer={}
    r_dum={}
    r_rand={}
    num_=[]
    infer=[]
    dum=[]
    rand=[]
    label=""
    for line in file: 
        if line.startswith("Sample"):
            if len(infer)!=0:
                r_infer[label]=infer
                r_dum[label]=dum
                r_rand[label]=rand
                num_infer[label]=num_
                infer=[]
                dum=[]
                rand=[]
                num_=[]
            label=line.rstrip()
        else:
            info=line.split(",")
            num_.append(int(info[0]))
            infer.append(float(info[1]))
            dum.append(float(info[2]))
            rand.append(float(info[3]))
    
    r_infer[label]=infer
    r_dum[label]=dum
    r_rand[label]=rand
    num_infer[label]=num_
    for label in r_infer.keys():
        print(label)
        avg_num=sum(num_infer[label])/len(num_infer[label])
        avg_r_inf=sum(r_infer[label])/len(r_infer[label])
        avg_r_dum=sum(r_dum[label])/len(r_dum[label])
        avg_r_rand=sum(r_rand[label])/len(r_rand[label])
        print("Avg number of cells: %2.1f"%avg_num)
        print("Avg R value from LP inference: %5.4f"%avg_r_inf)
        print("Avg R value by dummy: %5.4f"%avg_r_dum)
        print("Avg R value by random synthetic cell lines: %5.4f\n"%avg_r_rand)


                
            

if __name__ == "__main__":
    main()