import torch.nn
from param import *


class C3F2(torch.nn.Module):
    def __init__(self, in_sz=510, in_chl=2, chl_inc=4, c1k_sz=7, c2k_sz=7, c3k_sz=7, pk_sz=7, pst=1, fc1_out=160,
                 out_sz=10, verbose=False):
        super(C3F2, self).__init__()
        pk1_sz = pk_sz
        c1o_chl = in_chl + chl_inc
        c2o_chl = c1o_chl + chl_inc
        c3o_chl = c2o_chl + chl_inc
        c1o_sz = in_sz - c1k_sz + 1
        p1o_sz = int((c1o_sz - pk1_sz) / pst + 1)
        c2o_sz = p1o_sz - c2k_sz + 1
        p2o_sz = int((c2o_sz - pk_sz) / pst + 1)
        c3o_sz = p2o_sz - c3k_sz + 1
        fc1_in = c3o_chl * c3o_sz

        self.conv1 = torch.nn.Conv1d(in_channels=in_chl, out_channels=c1o_chl, kernel_size=c1k_sz)
        self.pool1 = torch.nn.AvgPool1d(kernel_size=pk1_sz, stride=pst)
        self.conv2 = torch.nn.Conv1d(in_channels=c1o_chl, out_channels=c2o_chl, kernel_size=c2k_sz)
        self.pool2 = torch.nn.AvgPool1d(kernel_size=pk_sz, stride=pst)
        self.conv3 = torch.nn.Conv1d(in_channels=c2o_chl, out_channels=c3o_chl, kernel_size=c3k_sz)
        self.fc1 = torch.nn.Linear(in_features=fc1_in, out_features=fc1_out)
        self.fc2 = torch.nn.Linear(in_features=fc1_out, out_features=out_sz)
        if verbose:
            print(f"c1_size: {in_sz},  ck1_sz: {c1k_sz},  in_chl: {in_chl}")
            print(f"c1o_sz: {c1o_sz},  c1o_chl: {c1o_chl}")
            print(f"pk_sz: {pk_sz},  p1o_size: {p1o_sz}")
            print(f"c2o_sz: {c2o_sz},  c2o_chl: {c2o_chl}")
            print(f"fc1_in: {fc1_in},  fc1_out: {fc1_out},  out_size: {out_sz}", flush=True)

    def forward(self, x):
        x = self.conv1(x)
        x = torch.relu(x)
        x = self.pool1(x)
        x = torch.relu(self.conv2(x))
        x = self.pool2(x)
        x = torch.relu(self.conv3(x))
        x = x.view(x.shape[0], x.shape[1] * x.shape[2])
        x = torch.relu(self.fc1(x))
        x = self.fc2(x)
        return x
