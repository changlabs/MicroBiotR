# Key feature diagram sources

The README figure uses the fixed-size vector layout in [key-features.svg](../figure/key-features.svg), exported to PNG for display. The Mermaid definitions below describe the workflow connections; their automatic layout may differ from the figure.

Blue parallelograms show user inputs, amber ovals show functions or scripts, and green rectangles show outputs. Each diagram follows a short workflow from left to right. An output named in one diagram can be supplied as input to another.

### Create shared flow-cytometry features

Start with raw FCS files and gating, or supply the downloaded gated events directly to `MBR_som`.

```mermaid
flowchart LR
    raw[/"<b>USER INPUT</b><br/>Raw FCS files"/] --> gating(["<b>FUNCTION</b><br/>Gating workflow"])
    gating --> events["<b>OUTPUT</b><br/>Gated event list"]
    events --> som(["<b>FUNCTION</b><br/>Step 1: MBR_som"])
    gated[/"<b>USER INPUT</b><br/>Downloaded gated events"/] --> som
    som --> abundance["<b>OUTPUT</b><br/>Cluster abundance table"]
    som --> codebook["<b>OUTPUT</b><br/>SOM prototype table"]
    som --> mapped["<b>OUTPUT</b><br/>Mapped FCS files"]
    classDef input fill:#DBEAFE,stroke:#1D4ED8,color:#172554,stroke-width:2px;
    classDef function fill:#FEF3C7,stroke:#B45309,color:#451A03,stroke-width:2px;
    classDef output fill:#DCFCE7,stroke:#15803D,color:#052E16,stroke-width:2px;
    class raw,gated input;
    class gating,som function;
    class events,abundance,codebook,mapped output;
```

### Compare groups and select features

Use the SOM abundance table or a 16S/RNA feature table, together with matched sample metadata.

```mermaid
flowchart LR
    statsin[/"<b>USER INPUT</b><br/>Feature table + metadata<br/>(steps 1, 14 or 15)"/] --> stat(["<b>FUNCTION</b><br/>Step 2: MBR_stat"])
    stat --> tested["<b>OUTPUT</b><br/>P-values and<br/>significant features"]
    tested --> groupplot(["<b>FUNCTION</b><br/>Steps 3–4: MBR_circle / MBR_violin"])
    groupplot --> groupfig["<b>OUTPUT</b><br/>Group abundance plots"]
    betain[/"<b>USER INPUT</b><br/>Feature table + metadata<br/>(steps 1, 14 or 15)"/] --> beta(["<b>FUNCTION</b><br/>Step 5: MBR_beta"])
    beta --> diversity["<b>OUTPUT</b><br/>PCoA and PERMANOVA"]
    fsin[/"<b>USER INPUT</b><br/>Feature table + metadata<br/>(steps 1, 14 or 15)"/] --> fs(["<b>FUNCTION</b><br/>Step 6: MBR_fs"])
    fs --> selected["<b>OUTPUT</b><br/>Selected features<br/>and importance plots"]
    classDef input fill:#DBEAFE,stroke:#1D4ED8,color:#172554,stroke-width:2px;
    classDef function fill:#FEF3C7,stroke:#B45309,color:#451A03,stroke-width:2px;
    classDef output fill:#DCFCE7,stroke:#15803D,color:#052E16,stroke-width:2px;
    class statsin,betain,fsin input;
    class stat,groupplot,beta,fs function;
    class tested,groupfig,diversity,selected output;
```

### Explore selected features and prototype profiles

Supply the selected features with the metadata or SOM prototype information needed by each analysis.

```mermaid
flowchart LR
    heatin[/"<b>USER INPUT</b><br/>Selected features<br/>+ SOM prototypes"/] --> heatmap(["<b>FUNCTION</b><br/>Step 7: MBR_heatmap"])
    heatmap --> heatout["<b>OUTPUT</b><br/>Channel-profile heatmap"]
    mantelin[/"<b>USER INPUT</b><br/>Selected features<br/>+ numeric metadata"/] --> mantel(["<b>FUNCTION</b><br/>Step 8: MBR_mantel"])
    mantel --> mantelout["<b>OUTPUT</b><br/>Correlation and Mantel plot"]
    modelin[/"<b>USER INPUT</b><br/>Selected features<br/>+ group metadata"/] --> model(["<b>FUNCTION</b><br/>Steps 9–10: MBR_ml / MBR_conf"])
    model --> modelout["<b>OUTPUT</b><br/>ROC and confusion plots"]
    clusterin[/"<b>USER INPUT</b><br/>SOM prototypes<br/>+ number of clusters"/] --> recluster(["<b>FUNCTION</b><br/>Step 11: MBR_reclustering"])
    recluster --> clusterout["<b>OUTPUT</b><br/>Prototype cluster labels"]
    classDef input fill:#DBEAFE,stroke:#1D4ED8,color:#172554,stroke-width:2px;
    classDef function fill:#FEF3C7,stroke:#B45309,color:#451A03,stroke-width:2px;
    classDef output fill:#DCFCE7,stroke:#15803D,color:#052E16,stroke-width:2px;
    class heatin,mantelin,modelin,clusterin input;
    class heatmap,mantel,model,recluster function;
    class heatout,mantelout,modelout,clusterout output;
```

### Inspect and export flow populations

Read a mapped FCS sample and prepare the selected prototypes for channel plots. Export selected events with matching FCS metadata.

```mermaid
flowchart LR
    fcsin[/"<b>USER INPUT</b><br/>Mapped FCS files<br/>+ selected sample"/] --> read(["<b>FUNCTION</b><br/>Step 12: MBR_read / MBR_process"])
    read --> events["<b>OUTPUT</b><br/>Processed event table"]
    protoin[/"<b>USER INPUT</b><br/>SOM prototypes<br/>+ selected cluster IDs"/] --> prepare(["<b>FUNCTION</b><br/>Step 12: MBR_prepare"])
    prepare --> bins["<b>OUTPUT</b><br/>Selected prototype table"]
    events --> plot(["<b>FUNCTION</b><br/>Step 12: MBR_flow_plot / MBR_plot"])
    bins --> plot
    plot --> plotout["<b>OUTPUT</b><br/>Channel plots"]
    exportin[/"<b>USER INPUT</b><br/>Selected events<br/>+ matching FCS metadata"/] --> save(["<b>FUNCTION</b><br/>Step 13: MBR_save"])
    save --> exportout["<b>OUTPUT</b><br/>Filtered FCS file"]
    classDef input fill:#DBEAFE,stroke:#1D4ED8,color:#172554,stroke-width:2px;
    classDef function fill:#FEF3C7,stroke:#B45309,color:#451A03,stroke-width:2px;
    classDef output fill:#DCFCE7,stroke:#15803D,color:#052E16,stroke-width:2px;
    class fcsin,protoin,exportin input;
    class read,prepare,plot,save function;
    class events,bins,plotout,exportout output;
```

