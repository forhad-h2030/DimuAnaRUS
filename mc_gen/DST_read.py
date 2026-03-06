#!/usr/bin/env python3
"""
Read E1039 DST ROOT files using PyROOT with shared libraries.
This mimics the C++ ReadDstTree function.
"""

import sys
import ROOT

def ReadDstTree(fn_dst, n_evt_max=0):
    """
    Read DST tree from E1039 ROOT file.
    
    Parameters:
    -----------
    fn_dst : str
        Path to DST ROOT file
    n_evt_max : int
        Maximum number of events to process (0 = all)
    """
    
    # Load required libraries
    print("Loading shared libraries...")
    try:
        ROOT.gSystem.Load('libinterface_main.so')
        print("  libinterface_main.so loaded")
    except Exception as e:
        print(f"Error loading libinterface_main.so: {e}")
        print("\nMake sure the library is in your LD_LIBRARY_PATH")
        print("You may need to set it up with:")
        print("  export LD_LIBRARY_PATH=/path/to/libraries:$LD_LIBRARY_PATH")
        sys.exit(1)
    
    # Try loading additional libraries if needed
    try:
        ROOT.gSystem.Load('libktracker.so')
        print("  libktracker.so loaded")
    except:
        pass  # Optional library
    
    # Open the ROOT file
    print(f"\nOpening file: {fn_dst}")
    file = ROOT.TFile.Open(fn_dst)
    if not file or file.IsZombie():
        print(f"Cannot open '{fn_dst}'. Abort.")
        sys.exit(1)
    
    print(f"File opened successfully: {file.GetName()}")
    
    # =========================================================================
    # Get the RUN node
    # =========================================================================
    print("\n" + "="*70)
    print("READING RUN TREE")
    print("="*70)
    
    tree_run = file.Get("T1")
    if not tree_run:
        print("Cannot get the RUN tree. Abort.")
        sys.exit(1)
    
    n_run_entries = tree_run.GetEntries()
    print(f"RUN tree entries: {n_run_entries}")
    
    if n_run_entries != 1:
        print(f"Find {n_run_entries} entries (!= 1) in the RUN tree. Abort.")
        sys.exit(1)
    
    # Set branch address for SQRun
    sq_run = ROOT.SQRun()
    tree_run.SetBranchAddress("RUN.SQRun", sq_run)
    tree_run.GetEntry(0)
    
    run_id = sq_run.get_run_id()
    print(f"Run ID = {run_id}")
    
    # =========================================================================
    # Get the DST node
    # =========================================================================
    print("\n" + "="*70)
    print("READING DST TREE")
    print("="*70)
    
    tree = file.Get("T")
    if not tree:
        print("Cannot get the DST tree. Abort.")
        sys.exit(1)
    
    n_evt = tree.GetEntries()
    print(f"N of events = {n_evt}")
    
    if n_evt_max != 0 and n_evt > n_evt_max:
        n_evt = n_evt_max
        print(f"Processing only first {n_evt} events")
    
    # Set branch addresses
    sq_evt = ROOT.SQEvent()
    sq_hit_vec = ROOT.SQHitVector()
    tree.SetBranchAddress("DST.SQEvent", sq_evt)
    tree.SetBranchAddress("DST.SQHitVector", sq_hit_vec)
    
    print("\n" + "="*70)
    print("PROCESSING EVENTS")
    print("="*70 + "\n")
    
    # Loop over events
    for i_evt in range(n_evt):
        tree.GetEntry(i_evt)
        
        # Print class names for first event
        if i_evt == 0:
            print(f"Class {sq_evt.ClassName()} {sq_hit_vec.ClassName()}\n")
        
        # Get event information
        evt_id = sq_evt.get_event_id()
        trig_bits = sq_evt.get_trigger()
        n_hit = sq_hit_vec.size()
        
        # Print event info
        print(f"E {i_evt} {evt_id} {trig_bits} {n_hit}")
        
        # Loop over hits
        for i_hit in range(n_hit):
            hit = sq_hit_vec.at(i_hit)
            hit_id = hit.get_hit_id()
            det_id = hit.get_detector_id()
            ele_id = hit.get_element_id()
            tdc = hit.get_tdc_time()
            
            # Print hit info
            print(f"  H {hit_id} {det_id} {ele_id} {tdc:.2f}")
        
        # Stop after 10 events (matching C++ code)
        if i_evt >= 10:
            print("\nStopping after 10 events (as in original C++ code)")
            break
    
    print("\n" + "="*70)
    print("DONE")
    print("="*70)
    
    file.Close()


def ReadDstTreeDetailed(fn_dst, n_evt_max=0, max_hits_per_event=None):
    """
    Read DST tree with more detailed output.
    
    Parameters:
    -----------
    fn_dst : str
        Path to DST ROOT file
    n_evt_max : int
        Maximum number of events to process (0 = all)
    max_hits_per_event : int
        Maximum number of hits to print per event (None = all)
    """
    
    # Load libraries
    print("Loading shared libraries...")
    ROOT.gSystem.Load('libinterface_main.so')
    
    # Open file
    print(f"\nOpening file: {fn_dst}")
    file = ROOT.TFile.Open(fn_dst)
    if not file or file.IsZombie():
        print(f"Cannot open '{fn_dst}'. Abort.")
        sys.exit(1)
    
    # Get RUN info
    tree_run = file.Get("T1")
    if tree_run:
        sq_run = ROOT.SQRun()
        tree_run.SetBranchAddress("RUN.SQRun", sq_run)
        tree_run.GetEntry(0)
        run_id = sq_run.get_run_id()
        print(f"\nRun ID: {run_id}")
    
    # Get DST tree
    tree = file.Get("T")
    if not tree:
        print("Cannot get the DST tree. Abort.")
        sys.exit(1)
    
    n_evt = tree.GetEntries()
    print(f"Total events: {n_evt}\n")
    
    if n_evt_max != 0 and n_evt > n_evt_max:
        n_evt = n_evt_max
    
    # Set branch addresses
    sq_evt = ROOT.SQEvent()
    sq_hit_vec = ROOT.SQHitVector()
    tree.SetBranchAddress("DST.SQEvent", sq_evt)
    tree.SetBranchAddress("DST.SQHitVector", sq_hit_vec)
    
    # Statistics
    total_hits = 0
    
    # Process events
    for i_evt in range(n_evt):
        tree.GetEntry(i_evt)
        
        evt_id = sq_evt.get_event_id()
        spill_id = sq_evt.get_spill_id()
        run_id = sq_evt.get_run_id()
        trig_bits = sq_evt.get_trigger()
        n_hit = sq_hit_vec.size()
        total_hits += n_hit
        
        print(f"\n{'='*70}")
        print(f"Event {i_evt}")
        print(f"{'='*70}")
        print(f"  Run ID:    {run_id}")
        print(f"  Event ID:  {evt_id}")
        print(f"  Spill ID:  {spill_id}")
        print(f"  Trigger:   {trig_bits}")
        print(f"  N hits:    {n_hit}")
        
        if n_hit > 0:
            print(f"\n  {'Hit':<6} {'TrackID':<10} {'DetID':<8} {'ElemID':<8} {'TDC Time':<12} {'Drift Dist':<12}")
            print(f"  {'-'*66}")
            
            n_print = n_hit if max_hits_per_event is None else min(n_hit, max_hits_per_event)
            
            for i_hit in range(n_print):
                hit = sq_hit_vec.at(i_hit)
                hit_id = hit.get_hit_id()
                track_id = hit.get_track_id()
                det_id = hit.get_detector_id()
                ele_id = hit.get_element_id()
                tdc = hit.get_tdc_time()
                drift = hit.get_drift_distance()
                
                print(f"  {hit_id:<6} {track_id:<10} {det_id:<8} {ele_id:<8} {tdc:<12.2f} {drift:<12.4f}")
            
            if n_hit > n_print:
                print(f"  ... and {n_hit - n_print} more hits")
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"Total events processed: {n_evt}")
    print(f"Total hits: {total_hits}")
    print(f"Average hits per event: {total_hits/n_evt:.1f}")
    print(f"{'='*70}\n")
    
    file.Close()


# =========================================================================
# Main execution
# =========================================================================
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Read E1039 DST ROOT files using PyROOT',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s /data2/e1039/dst/run_004621_spin.root
  %(prog)s /data2/e1039/dst/run_004621_spin.root -n 50
  %(prog)s /data2/e1039/dst/run_004621_spin.root --detailed
  %(prog)s /data2/e1039/dst/run_004621_spin.root --detailed --max-hits 20

Note: This script requires libinterface_main.so to be in your LD_LIBRARY_PATH
        """
    )
    
    parser.add_argument('filename', 
                       help='Path to DST ROOT file')
    parser.add_argument('-n', '--nevents', type=int, default=0,
                       help='Maximum number of events to process (0=all, default=0)')
    parser.add_argument('--detailed', action='store_true',
                       help='Use detailed output format')
    parser.add_argument('--max-hits', type=int, default=None,
                       help='Maximum hits to display per event (for detailed mode)')
    
    args = parser.parse_args()
    
    try:
        if args.detailed:
            ReadDstTreeDetailed(args.filename, args.nevents, args.max_hits)
        else:
            ReadDstTree(args.filename, args.nevents)
    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        sys.exit(0)
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
