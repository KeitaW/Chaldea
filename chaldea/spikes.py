# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 17:09:37 2014

@author: keitawatanabe
"""

from brian import *
import json
import matplotlib.pyplot as plt
import pylab
import pickle
import os
from itertools import cycle
class Spikes:
    '''
    人工的なスパイクトレインを生成するクラス．init時にnumOfNeuron個のニューロンが周波数frequencyでポアソン発火した際のスパイクトレインを生成する.
    なお，スパイクトレインはms単位で測ったものである．
    '''
    def __init__(self, numOfNeurons, frequency, runTime):
        self.num_of_neurons = numOfNeurons
        self.frequency = frequency
        self.run_time = runTime
        print( "run_time",self.run_time)
        # シミュレーションの時間解像度，clockを1msおきに変更．（実験データの方に合わせるため)
        self.unit_time = 1*ms
        my_clock=Clock(dt=self.unit_time)
        P = PoissonGroup(self.num_of_neurons, frequency, my_clock)
        M = SpikeMonitor(P)
        reinit_default_clock(t=0.0*second)  # simのクロックが必ず0から始まるようにする
        run(self.run_time)
        self.spikes = [[float(i[1]/second),i[0], 0] for i in M.spikes]      # スパイクトレインをリストの形で格納（第三引数は色の識別用（青がポアソンスパイク，赤が人工スパイク））
        print( 'finished spikes generation')

    def sequence_generator(reference_neuron, length_of_sequence,duration_of_sequence, num_of_sequence, total_num_of_neurons,  time, jitter=0.02):
        '''
        generate random sequence
        arguments:
        '''

        appear_times = np.sort(np.random.randint(low=0, high=time, size=num_of_sequence))
        sequence_neurons = np.insert(np.random.randint(low=0, high=total_num_of_neurons, size=length_of_sequence-1), 0, reference_neuron)
        sequence_times = np.sort(np.random.random_sample(size=length_of_sequence)) * duration_of_sequence
        sequence = []
        for atime in appear_times:
            # todo: jitterをはしご式に伸ばしていけるように改良すること
            jitters = np.random.random_sample(size=length_of_sequence) * jitter
            t = atime + sequence_times + jitters
            t[0] = atime        # reference neuronがsequenceの先頭に来るための処理
            sequence_seed = [[i, j] for i,j in zip(t[1:], sequence_neurons)]
            sequence += sequence_seed

        return sequence, sequence_neurons

    def add_spikes(self, added_spikes, color_id = 1):
        '''
        listで与えられるようなスパイクトレインをself.spikesに加える．listは（spikeTime, neuron）の形で渡すこと．
        なお，スパイクトレインは発火時間でソートされる．
        '''
        self.added_spikes = added_spikes
        self.spikes += [list([i[0], i[1], color_id]) for i in added_spikes]
        self.spikes.sort(key=lambda x:(x[0]))

    def extract_and_padding(self, raw_fragments, dt):
        '''
        下記create_fragments関数にて生成されたraw_fragments関数からfragments文字列
        '''
        fragments = []
        starting_point_of_fragments = []
        for fragment in raw_fragments:
            start_time = fragment[0,0]
            sequence_id = fragment[0,2]
            starting_point_of_fragments.append([start_time, sequence_id])
            end_time = fragment[-1,0]
            times = np.arange(start_time, start_time+0.5, dt)
            print('len_times', len(times))
            candidate_contanor = np.array([],dtype='str')
            for t in times:
                candidate_idx = np.where(np.logical_and(np.abs(fragment[:,0]-t)<=dt, fragment[:,0]-t>=0))[0]
                candidate = fragment[candidate_idx, 1].astype(int)
                candidate_str = candidate if len(candidate)!=0 else '-'
                candidate_contanor = np.append(candidate_contanor, candidate_str)
            fragments.append(candidate_contanor)
        return fragments, starting_point_of_fragments


    def extract(self, raw_fragments, dt):
        '''
        下記create_fragments関数にて生成されたraw_fragments関数からfragments文字列を生成padding無し
        '''
        fragments = []
        starting_point_of_fragments = []
        for fragment in raw_fragments:
            start_time = fragment[0,0]
            starting_point_of_fragments.append(start_time)
            end_time = fragment[-1,0]
            times = np.arange(start_time, end_time, dt)
            candidate_contanor = np.array([],dtype='str')
            for t in times:
                candidate_idx = np.where(np.logical_and(np.abs(fragment[:,0]-t)<=dt, fragment[:,0]-t>=0))[0]
                candidate = fragment[candidate_idx, 1].astype(int)
                if len(candidate)!=0:
                    candidate_contanor = np.append(candidate_contanor, candidate)
            fragments.append(candidate_contanor)
        return fragments

    def create_fragments(self, referenceneuron, exduration=0.5, dt=0.001):
        '''
        reference neuronにとったニューロンの各発火から，timelen*dt秒分の発火を抽出した後，それらをreference neuronの発火で分割．それらを起点にして他のニューロンの発火をidのみを記録する形で表現する．
        exdurationは各ビンの大きさを規定している．
        '''
        spikes = np.array(self.spikes)
        # 各reference neuronの発火時間を示すindexを抽出
        spoints = np.where(spikes[:,1]==referenceneuron)[0] # 正しい
        stimes = [ spike[0] for spike in spikes[spoints]] # 正しい．がリスト内包表記だけで済むようにもおもふ
        extract_indices = [np.where(np.logical_and(spikes[:,0]>=time, spikes[:,0]<=time+exduration))[0] for time in stimes]
        # 空の要素を削除 -> 同じneuron id のものしか抽出されていない．要改善
        length_of_shortest_fragment = 3
        indices = [eindices for eindices in extract_indices if eindices.size >= length_of_shortest_fragment]
        #self.fragments = [spikes[idx,1] for idx in indices]
        self.raw_fragments = [spikes[idx] for idx in indices]
        self.fragments, self.starting_point_of_fragments = self.extract_and_padding(self.raw_fragments,dt)
        self.fragments_without_padding = self.extract(self.raw_fragments, dt)

    def save_fragments(self, dir_str):
        '''
        self.fragmentsを保存するための関数．本来はsave_resultにまとめるべきかとも思われるが，今回は後付でスパイクシーケンスを導入したりするのでこの形にした．
        '''
        wf = open(dir_str + '/raw_fragments.npy', 'w')
        np.save(wf, self.raw_fragments)
        wf.close()
        wf = open(dir_str + '/raw_fragments.ssv', 'w')
        np.savetxt(wf, self.raw_fragments, fmt="%s")
        wf.close()

        wf = open(dir_str + '/raw_fragments_without_padding.npy', 'w')
        np.save(wf, self.fragments_without_padding)
        wf.close()
        wf = open(dir_str + '/raw_fragments_without_padding.ssv', 'w')
        np.savetxt(wf, self.fragments_without_padding, fmt="%s")
        wf.close()


        wf = open(dir_str + '/fragments.ssv', 'w')
        np.savetxt(wf, self.fragments, fmt="%s")
        wf.close()
        wf = open(dir_str + '/fragments.npy', 'w')
        np.save(wf, self.fragments)
        wf.close()
        # 以下starting_point_of_fragmentsの保存
        wf = open(dir_str + '/starting_points_of_fragments.npy', 'w')
        np.save(wf, self.starting_point_of_fragments)
        wf.close()

    def save_result(self, dir_str):
        '''
        シミュレーションの結果をdir_strで指定されたディレクトリに保存する．
        '''
        # スパイクシミュレーションデータの保存
        wf = open(dir_str + "/spikes.ssv", 'w')
        np.savetxt(wf,self.spikes)     # scientific reportのデータと書式を合わせている．
        wf.close()
        # スパイクシミュレーションデータの保存
        wf = open(dir_str + "/spikes.npy", 'w')
        np.save(wf,self.spikes)     # scientific reportのデータと書式を合わせている．
        wf.close()

        wf = open(dir_str + "/summery.json", 'w')
        summery = {"number of neurons":self.num_of_neurons,  "frequency":self.frequency, "runTime":self.run_time}
        json.dump(summery, wf , indent=4)
        wf.close()
        # 以下スパイクデータのプロット
        plt.title("raster plot of spikes")
        plt.xlabel("time (s)")
        plt.ylabel("neuron")
        colors = 'crbgmky'
        length = int(self.run_time)
        # plt.figure(figsize=(length, 4))
        for spike in self.spikes:
            plt.plot(spike[0], spike[1], colors[spike[2]] + '.')
        wf = open(dir_str + "/spikes.png", 'w')
        plt.savefig(wf, format='png')
        plt.clf()
        wf.close()

        # 以下このクラスそのものをpickleで保存
        wf = open(dir_str + "/spikes_class.dat", 'w')
        pickle.dump(self, wf)
        wf.close()
