#!/usr/bin/env python3
import os, json, argparse, sys, datetime
from pyliftover import LiftOver

URN  = os.environ.get("MAVEDB_URN", "na")
STEP = os.environ.get("STEP", "liftover")
def log(reason, subid="na", **kv):
    ts = datetime.datetime.now().astimezone().isoformat()
    extras = " ".join(f"{k}={json.dumps(v)}" for k, v in kv.items())
    sys.stderr.write(f"[{ts}][MaveDB][URN={URN}][STEP={STEP}][REASON={reason}][SUBID={subid}] {extras}\n")
    sys.stderr.flush()

def main():
  parser = argparse.ArgumentParser(description='Lift-over variants associated with MaveDB scores')
  parser.add_argument('--metadata', type=str, help="path to file with MaveDB metadata")
  parser.add_argument('--mapped_variants', type=str, help="path to file with variants mapped to MaveDB scores")
  parser.add_argument('--reference', type=str, default="hg38", help="genome (default: 'hg38')")
  args = parser.parse_args()

  reference = args.reference
  mapped    = args.mapped_variants

  # Load metadata
  metadata  = json.load(open(args.metadata))

  genome = metadata.get('extraMetadata', {}).get('reference', 'hg38')
  log("liftover_start", genome=genome, reference=reference, mapped=mapped)

  if genome == reference:
    # just rename file if variants are already mapped to reference genome
    new_name = f"liftover_{mapped}"
    os.rename(mapped, new_name)
    log("rename_passthrough", out=new_name)
  else:
    liftover_variants(mapped, genome, reference)

  return True

def liftover_variants(mapped, genome, reference):
  # write file information with lifted-over coordinates to new file
  log("converting", genome=genome, reference=reference)
  chain = LiftOver(f"{genome}To{reference.capitalize()}.over.chain.gz")

  out_path = f"liftover_{mapped}"
  out = open(out_path, "w")

  n_total = n_lifted = 0

  with open(mapped) as f:
    header = f.readline()
    out.write(header)

    # Try to infer the accession column for SUBID; fallback to coord token
    header_cols = [c.strip() for c in header.rstrip("\n").split("\t")]
    try:
      acc_idx = header_cols.index("accession")
    except ValueError:
      acc_idx = None

    for line in f:
      n_total += 1
      l = line.rstrip("\n").split("\t")
      chr_, start, end = l[0:3]

      # preferred SUBID: accession; else coord token
      subid = l[acc_idx] if acc_idx is not None and acc_idx < len(l) else f"{chr_}:{start}-{end}"

      # convert start
      conv = chain.convert_coordinate("chr" + chr_, int(start))
      if len(conv) == 0:
        log("liftover_unmapped_start", subid=subid, chr=chr_, pos=start)
        raise Exception("unmapped coordinate")
      if len(conv) > 1:
        log("liftover_multiple_mappings_start", subid=subid, chr=chr_, pos=start, n=len(conv))
        raise Exception("multiple coordinates returned")
      new_chr, new_start, strand, size = conv[0]

      # convert end
      conv = chain.convert_coordinate("chr" + chr_, int(end))
      if len(conv) == 0:
        log("liftover_unmapped_end", subid=subid, chr=chr_, pos=end)
        raise Exception("unmapped coordinate")
      if len(conv) > 1:
        log("liftover_multiple_mappings_end", subid=subid, chr=chr_, pos=end, n=len(conv))
        raise Exception("multiple coordinates returned")
      new_chr, new_end, strand, size = conv[0]

      # write lifted row
      l[0:3] = new_chr.replace("chr", ""), str(new_start), str(new_end)
      out.write('\t'.join(l) + "\n")
      n_lifted += 1

  out.close()
  log("liftover_done", processed=n_total, lifted=n_lifted, out=out_path)

if __name__ == "__main__":
  main()