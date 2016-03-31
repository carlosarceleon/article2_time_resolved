def publish():
    import os
    import shutil

    article_folder = "/home/carlos/Documents/PhD/Articles/Article_2/4113877skpkyw/"
    article_files = [f for f in os.listdir(article_folder) if os.path.splitext(f)[1] == '.png']

    plot_directories = [
        #'./SpectrumResults/',
        './Results/',
        './Results_v2/'
    ]

    for pd in plot_directories:
        local_files = [
            f for f in os.listdir(pd) if os.path.splitext(f)[1] == '.png'
        ]

        for lf in local_files:
           if lf in article_files:
              print "   Copying {0}".format(lf)
              shutil.copy(os.path.join(pd,lf),article_folder)
