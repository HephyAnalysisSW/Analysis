import ROOT, os
import keras
import Analysis.Tools.syncer
import RootTools.plot.helpers as plot_helpers

def plot_metrics(history, plot_directory):

    # Make sure the plot dir exists and copy index.php
    if not os.path.exists(plot_directory):
        os.makedirs(plot_directory)
    plot_helpers.copyIndexPHP( plot_directory )

    # Get all keys (only store strings without 'val_')
    keys = []
    for key in list(history.history.keys()):
        if "val_" in key:
            continue
        else:
            keys.append(key)

    # Loop over keys
    for key in keys:
        # security check
        if 'val_'+key not in list(history.history.keys()):
            print('[ERROR] val_'+key+' not in history.keys(). The '+key+' metric will not be plotted.')
            continue

        # Get the numbers
        data_train = history.history[key]
        data_val = history.history['val_'+key]
        N_epochs_train = len(data_train)
        N_epochs_val = len(data_val)

        # create graphs for train and validation sets
        graph_train = ROOT.TGraph(N_epochs_train)
        graph_val = ROOT.TGraph(N_epochs_val)
        for i in range(N_epochs_train): graph_train.SetPoint(i+1, i+1, data_train[i])
        for i in range(N_epochs_val):   graph_val.SetPoint(i+1, i+1, data_val[i])

        # Plot
        c = ROOT.TCanvas(key, key, 600, 600)
        graph_train.SetLineWidth(3)
        graph_train.SetLineColor(ROOT.kRed-3)
        graph_val.SetLineWidth(3)
        graph_val.SetLineColor(ROOT.kAzure-7)
        graph_train.Draw("AC") # has to be here for TGraph to have an axis object

        # Set Labels, range and draw again
        graph_train.GetXaxis().SetTitle("Epoch")
        graph_train.GetYaxis().SetTitle(key)
        maximum = max(data_train) if max(data_train) > max(data_val) else max(data_val)
        graph_train.GetYaxis().SetRangeUser(0.,1.8*maximum)
        graph_train.Draw("AC")
        c.Update();
        graph_val.Draw("SAME C")

        # Legend
        leg = ROOT.TLegend(0.6, 0.65, 0.85, 0.9)
        leg.AddEntry(graph_train, "Training", "l")
        leg.AddEntry(graph_val, "Validation", "l")
        leg.SetBorderSize(0)
        leg.Draw()

        # Save plot
        c.Print(plot_directory+"/"+key+".png")
