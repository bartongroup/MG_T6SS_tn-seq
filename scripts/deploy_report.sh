rsync -rvm --include='*/' --include='doc/report.html' --include='doc/image/*' --include='tab/*' \
  --exclude='*' . cluster:/cluster/gjb_lab/mgierlinski/public_html/ts66_tn_seq
